import numpy as np
from scipy.optimize import minimize
import scipy as sp
import os,sys
import time
import re
from statistics import mean
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from sklearn import *
import optimization as opt
import openbabel
import itertools
import symmetry as symm
from socket import gethostname
import shlex
from scipy.optimize import fmin
from scipy import optimize
import shutil
import copy
import traceback
import apicall as call
from scipy.interpolate import interp1d
from PyAstronomy import pyasl

def CheckAllVdwTypesExist(poltype,keyfilename):
    temp=open(keyfilename,'r')
    results=temp.readlines()
    temp.close()
    vdwtypes=[]
    for line in results:
        if 'vdw ' in line and '#' not in line:
            linesplit=line.split()
            atomtype=int(linesplit[1])
            vdwtypes.append(atomtype)
    alltypes=[]
    for atom in poltype.rdkitmol.GetAtoms():
        atomidx=atom.GetIdx()+1
        symtype=poltype.idxtosymclass[atomidx]
        alltypes.append(symtype)
    for atomtype in alltypes:
        if atomtype not in vdwtypes:
            raise ValueError('Missing vdwtype in keyfile '+str(atomtype)) 

def best_fit_slope_and_intercept(xs,ys):
    xs=np.array(xs)
    ys=np.array(ys)
    m = (((mean(xs)*mean(ys)) - mean(xs*ys)) /
         ((mean(xs)*mean(xs)) - mean(xs*xs)))
    b = mean(ys) - m*mean(xs)
    return m, b

def squared_error(ys_orig,ys_line):
    ys_orig=np.array(ys_orig)
    ys_line=np.array(ys_line)
    return sum((ys_line - ys_orig) * (ys_line - ys_orig))

def coefficient_of_determination(ys_orig,ys_line):
    y_mean_line = [mean(ys_orig) for y in ys_orig]
    squared_error_regr = squared_error(ys_orig, ys_line)
    squared_error_y_mean = squared_error(ys_orig, y_mean_line)
    try:
        value=1 - (squared_error_regr/squared_error_y_mean)

    except:
        value=1 # incase division by 0

    return value 

def MeanError(data,pred):
    Sum=0
    for i in range(len(data)):
        true=data[i]
        predv=pred[i]
        diff=true-predv
        Sum+=diff
    Sum=Sum/len(data)
    return Sum
    


def writePRM(poltype,params,vdwtypes,idxtotype):
    vdwradius=None
    vdwdepth=None
    vdwred=None
    vdwtypesalreadyfit=[]
    vdwtypetonewline={}
    for i in range(len(params)):
        oFile = open("temp.key", 'w')
        rFile=open(poltype.key4fname,'r')
        vdwtype=vdwtypes[i]
        vartype=idxtotype[i]
        vdwradius=None
        vdwdepth=None
        vdwred=1
        if vartype=='rad':
            vdwradius=params[i]

        elif vartype=='depth':
            vdwdepth=params[i]
            if (i+1) in idxtotype.keys():
                nexttype=idxtotype[i+1]

                if nexttype!='red':
                    vdwred=1
            else:
                vdwred=1

        elif vartype=='red':
            vdwred=params[i]

        for line in rFile.readlines():
            linesplit=line.split()
            if 'vdw' in line and '#' not in line:
                if len(linesplit)!=5:
                    linesplit.append('1')
                line=' '.join(linesplit)+'\n' 
                if linesplit[1]==vdwtype:
                    if vdwradius!=None:
                        linesplit[2]=str(vdwradius)
                    if vdwdepth!=None:
                        linesplit[3]=str(vdwdepth)
                    if vdwred!=None:
                        linesplit[4]=str(vdwred)
                    newline=' '.join(linesplit)+'\n' 
                    vdwtypetonewline[vdwtype]=newline
                    oFile.write(newline)
                else:
                    notinline=True
                    for othertype in vdwtypesalreadyfit:
                        if linesplit[1]==othertype:
                            notinline=False
                            break
                    if notinline==True:
                        oFile.write(line)
                    else:
                        newline=vdwtypetonewline[othertype]
                        oFile.write(newline)

            else:
                oFile.write(line)

        vdwtypesalreadyfit.append(vdwtype)
        vdwradius=None
        vdwdepth=None
        vdwred=None
    oFile.flush()
    os.fsync(oFile.fileno())
    oFile.close()
    rFile.close()
    os.remove(poltype.key4fname)
    os.rename("temp.key",poltype.key4fname)
    CheckAllVdwTypesExist(poltype,poltype.key4fname)
    return

def readOneColumn(filename,columnnumber,prefix=None):
    temp=open(filename,'r')
    results=temp.readlines()
    temp.close()
    array=[]
    for line in results:
        if prefix!=None:
            if prefix not in line:
                continue
        linesplit=line.split()
        value=linesplit[columnnumber]
        array.append(value)
    return array
     
def NormalizeTarget(poltype,filename,aprefix=None):
    prefixarray=readOneColumn(filename,0)
    energyarray=readOneColumn(filename,1)
    energyarray=[float(i) for i in energyarray]
    prefixtoenergyarray={}
    for prefixidx in range(len(prefixarray)):
        prefix=prefixarray[prefixidx]
        if aprefix!=None:
            if aprefix not in prefix:
                continue
        energy=energyarray[prefixidx]
        prefixsplit=prefix.split('_')
        prefixsplit=prefixsplit[:-1] 
        prefix=''.join(prefixsplit)
        if prefix not in prefixtoenergyarray.keys():
            prefixtoenergyarray[prefix]=[]
        prefixtoenergyarray[prefix].append(energy)
    for prefix,energyarray in prefixtoenergyarray.items():
        normenergyarray=[i-min(energyarray) for i in energyarray]
        prefixtoenergyarray[prefix]=normenergyarray
    totalenergyarray=[]
    for prefix,energyarray in prefixtoenergyarray.items():
        for e in energyarray:
            totalenergyarray.append(e)
    
    return np.array(totalenergyarray)


def myFUNC(params,poltype,vdwtypes,idxtotype,count):
    params=list(params)
    writePRM(poltype,params,vdwtypes,idxtotype)
    target = NormalizeTarget(poltype,'QM_DATA')

    temp=open('QM_DATA','r')
    cmdarray=[] 
    filenamearray=[]
    for line in temp.readlines():
        xyzname=line.split()[0]
        filename=xyzname.replace('.xyz','.alz')
        cmdstr=poltype.analyzeexe+' '+xyzname+' '+'-k '+poltype.key4fname +' e '+'> '+filename
        cmdarray.append(cmdstr)
        filenamearray.append(filename)

    temp.close()
    for cmdidx in range(len(cmdarray)):
        cmd=cmdarray[cmdidx]
        filename=filenamearray[cmdidx]
        poltype.call_subsystem([cmd],True)
        temp={cmd:filename} 
        finishedjobs,errorjobs=poltype.WaitForTermination(temp,False)

    ReadAnalyzeEnergiesWriteOut(poltype,filenamearray)
    current = NormalizeTarget(poltype,'SP.dat')

    current,target,distarray=ScreenHighEnergyPoints(poltype,current,target)
    if count>1:
        weightlist=np.exp(-np.array(target)/poltype.boltzmantemp)

    else:
        weightlist=np.ones(len(target))
    value=weightlist*(current-target)

    return weightlist*(current-target)

def ScreenHighEnergyPoints(poltype,current,target,distarray=None):
    tol=10
    newcurrent=[]
    newtarget=[]
    newdistarray=[]
    for i in range(len(current)):
        cur=current[i]
        tar=target[i]
        try:
            if distarray==None:
                pass
        except:
            dist=distarray[i]
        if tar>=tol:
            pass
        else:
            newcurrent.append(cur)
            newtarget.append(tar)
            try:
                if distarray==None:
                    pass
            except:
                newdistarray.append(dist)
    return np.array(newcurrent),np.array(newtarget),np.array(newdistarray)


def PlotQMVsMMEnergy(poltype,vdwtypesarray,prefix,count,classkeytofitresults,allprefix=False):
    if allprefix==False:
        target= NormalizeTarget(poltype,'QM_DATA',prefix)
        current=NormalizeTarget(poltype,'SP.dat',prefix)

    else:
        target=[]
        current=[]
        for ls in prefix:
            target.extend(NormalizeTarget(poltype,'QM_DATA',ls))
            current.extend(NormalizeTarget(poltype,'SP.dat',ls))
    current,target,distarray=ScreenHighEnergyPoints(poltype,current,target)
    vdwtypes=[str(i) for i in vdwtypesarray]
    vdwtypestring=','.join(vdwtypes)
    MSE=MeanError(current,target)
    MAE=metrics.mean_absolute_error(current,target)
    m, b = best_fit_slope_and_intercept(current,target)
    regression_line = [(m*x)+b for x in current]
    weightlist=np.exp(-np.array(target)/poltype.boltzmantemp)

    energy_list=target-current
    if count>1:
        energy_list=np.multiply(weightlist,energy_list)
    shifted_target=np.add(1,target)
    shifted_current=np.add(1,current)
    relative_energy_list=shifted_target-shifted_current
    if count>1:
        relative_energy_list=np.multiply(weightlist,relative_energy_list)

    def FirstRMSD(c):
        return np.sqrt(np.mean(np.square(np.add(target-current,c))))
    result=fmin(FirstRMSD,.5)
    first_new_rmse=FirstRMSD(0)


    def RMSD(c):
        return np.sqrt(np.mean(np.square(np.add(energy_list,c))))
    result=fmin(RMSD,.5)
    new_rmse=RMSD(0)

    def RMSDRel(c):
        return np.sqrt(np.mean(np.square(np.add(np.divide(relative_energy_list,shifted_target),c))))
    resultRel=fmin(RMSDRel,.5)
    new_rmse_rel=RMSDRel(0)
    r_squared = coefficient_of_determination(current,target)
    fig = plt.figure()
    plt.plot(current,target,'.',label='R^2=%s MAE=%s RMSE=%s RelRMSE=%s MSE=%s'%(round(r_squared,2),round(MAE,2),round(new_rmse,2),round(new_rmse_rel,2),round(MSE,2)))
    plt.plot(current,regression_line,label='Linear Regression line')
    plt.ylabel('QM BSSE Corrected (kcal/mol)')
    plt.xlabel('AMOEBA (kcal/mol)')
    plt.legend(loc='best')
    plt.title('QM vs AMOEBA , '+vdwtypestring)
    if count>1:
        suffix='_boltzman.png'
        ostring='Boltzmann Fit'
    else:
        suffix='.png'
        ostring=''
    if allprefix==False:
        fig.savefig('QMvsAMOEBA-'+prefix+'_'+vdwtypestring+suffix)
    else:
        prefstring=','.join(prefix)
        fig.savefig('QMvsAMOEBA-'+prefstring+'_'+vdwtypestring+suffix)
        prefix=prefstring
    for vdwtypenum in vdwtypes:
        thestring='VDW '+str(vdwtypenum)+' RMSD(MM,QM) '+str(new_rmse)+' '+'RelativeRMSD(MM,QM) '+str(new_rmse_rel)+' '+ostring+'\n'
        poltype.WriteToLog(thestring)
        classkeytofitresults[vdwtypenum]=thestring

    rmsetol=1.9
    relrmsetol=.2
    goodfit=True
    if new_rmse>rmsetol and new_rmse_rel>relrmsetol:
        goodfit=False
        if allprefix==True and count>1:
            raise ValueError('RMSE is too high! RMSE='+str(new_rmse)+' tol='+str(rmsetol)+' '+'RelRMSE='+str(new_rmse_rel)+' tol='+str(relrmsetol)+' '+prefix)
        elif allprefix==True and count<=1:
            mes='RMSE is too high! RMSE='+str(new_rmse)+' tol='+str(rmsetol)+' '+'RelRMSE='+str(new_rmse_rel)+' tol='+str(relrmsetol)+' '+prefix
            poltype.WriteToLog(mes)

    return goodfit,classkeytofitresults


def VDWOptimizer(poltype,count,fitredboolarray):
    x0 = []
    curvdwtypes=readOneColumn("INITIAL.PRM", 1)
    vdwtypes=[]
    temp=open("INITIAL.PRM",'r')
    lines = temp.readlines()
    idxtotype={}
    count=0
    rmax=readOneColumn("INITIAL.PRM", 6)
    rmax=[float(i) for i in rmax]
    rmin=readOneColumn("INITIAL.PRM", 5)
    rmin=[float(i) for i in rmin]
    depthmax=readOneColumn("INITIAL.PRM", 8)
    depthmax=[float(i) for i in depthmax]
    depthmin=readOneColumn("INITIAL.PRM", 7)
    depthmin=[float(i) for i in depthmin]
    redmax=readOneColumn("INITIAL.PRM", 10)
    redmax=[float(i) for i in redmax]
    redmin=readOneColumn("INITIAL.PRM", 9)
    redmin=[float(i) for i in redmin]
    l1=list(zip(rmin, rmax))
    l2=list(zip(depthmin, depthmax))
    l3=list(zip(redmin, redmax))
    lower=[]
    upper=[]
    rmintovdwtypes=dict(zip(rmin,curvdwtypes))  
    vdwtypestoprmsremove={}
    for i in range(len(l1)):
        fitredbool=fitredboolarray[i]
        l1bound=l1[i]
        l2bound=l2[i]
        l3bound=l3[i]
        l1lower=l1bound[0] 
        l1upper=l1bound[1] 
        l2lower=l2bound[0] 
        l2upper=l2bound[1] 
        l3lower=l3bound[0] 
        l3upper=l3bound[1]
        
        vdwtype=rmintovdwtypes[l1lower]
        if vdwtype not in vdwtypestoprmsremove.keys():
            vdwtypestoprmsremove[vdwtype]=[]
        if l1lower!=l1upper: 
            lower.append(l1lower)
        else:
            vdwtypestoprmsremove[vdwtype].append('S')
        if l2lower!=l2upper:
            lower.append(l2lower)
        else:
            vdwtypestoprmsremove[vdwtype].append('T')

        if fitredbool==True:
            lower.append(l3lower)
        if l1lower!=l1upper:
            upper.append(l1upper)
        if l2lower!=l2upper:
            upper.append(l2upper)
        if fitredbool==True:
            upper.append(l3upper)
    MyBounds=[lower,upper]
    MyBounds=tuple(MyBounds)
    for lineidx in range(len(lines)):
        vdwtype=curvdwtypes[lineidx]
        fitredbool=fitredboolarray[lineidx]
        line=lines[lineidx]
        prmstoremove=vdwtypestoprmsremove[vdwtype]
        radius=float(line.split()[2])
        depth=float(line.split()[3])
        red=float(line.split()[4])
        if 'S' not in prmstoremove:
            x0.append(radius)
            idxtotype[count]='rad'
            count+=1
            vdwtypes.append(vdwtype)

        if 'T' not in prmstoremove:
            x0.append(depth)
            idxtotype[count]='depth'
            count+=1
            vdwtypes.append(vdwtype)
        if fitredbool==True:
            x0.append(red)
            idxtotype[count]='red'
            count+=1
            vdwtypes.append(vdwtype)

    x0 = np.array(x0)
    ''' local optimization method can be BFGS, CG, Newton-CG, L-BFGS-B,etc.., see here\
    https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.optimize.minimize.html'''
    errorfunc= lambda p, poltype, vdwtypes,idxtotype,count: (myFUNC(p,poltype,vdwtypes,idxtotype,count))
    ret = optimize.least_squares(errorfunc, x0, jac='2-point', bounds=MyBounds, diff_step=1e-4,args=(poltype,vdwtypes,idxtotype,count))
    vdwradii=[]
    vdwdepths=[] 
    vdwreds=[]
    ofile = open("RESULTS.OPT", "a")
    for i in range(len(ret.x)):
        vartype=idxtotype[i]
        if vartype=='rad':
            vdwradius=round(ret.x[i],3)
            vdwradii.append(vdwradius)

        elif vartype=='depth':
            vdwdepth=round(ret.x[i],3)
            vdwdepths.append(vdwdepth)
            if (i+1) in idxtotype.keys():
                nexttype=idxtotype[i+1]
                if nexttype!='red':
                    vdwred=1
                    vdwreds.append(vdwred)
            else:
                vdwred=1
                vdwreds.append(vdwred)


        elif vartype=='red':
            vdwred=round(ret.x[i],3)
            vdwreds.append(vdwred)
            
    ofile.write("%s\n"%(ret.message))
    ofile.write("\n")
    ofile.close()
    return vdwradii,vdwdepths,vdwreds


def MoveDimerAboutMinima(poltype,txyzFile,outputprefixname,nAtomsFirst,atomidx1,atomidx2,equildistance,array):
    fnamearray=[]
    for i in range(len(array)):
        frac=array[i]
        fname=outputprefixname+'_%s.xyz' %(str(frac))
        fnamearray.append(fname)
        with open(fname,'w') as f:
            twoDots = [atomidx1,atomidx2] #two atoms selected to vary; the atom number exactly in TINKER file
            temp=open(txyzFile,'r')
            lines=temp.readlines()
            temp.close()
            coordSecondMol = [] #Coordinates of the second molecule
            data = lines[twoDots[0]].split()
            coordFirstAtm = [float(data[2]), float(data[3]), float(data[4])] 
            data = lines[twoDots[1]].split()
            coordSecondAtm = [float(data[2]), float(data[3]), float(data[4])] 
            varyVector = []
            norm = 0.0
            for i in range(3):
                varyVector.append(coordSecondAtm[i] - coordFirstAtm[i])
                norm += (coordSecondAtm[i]-coordFirstAtm[i])**2.0
            norm = norm**0.5 # initial distance between the two atoms
            distance=norm*frac
            for i in range(3):
                varyVector[i] = (distance - norm) * varyVector[i]/norm # scaling the x,y,z coordinates of displacement vector to new distance
            stringsList = [] 
            for i in range(nAtomsFirst+1, len(lines), 1):
                data = lines[i].split()
                coordSecondMol.append( [float(data[2]), float(data[3]), float(data[4]),] )
                headString = ' '.join(lines[i].split()[0:2]) 
                tailString = ' '.join(lines[i].split()[5:])
                stringsList.append([headString, tailString])
            for i in range(0, nAtomsFirst+1, 1):
                f.write(lines[i])
    
            for i in range(len(coordSecondMol)):
                for j in range(3):
                    coordSecondMol[i][j] = coordSecondMol[i][j]+varyVector[j]
    
            for i in range(len(coordSecondMol)):
                f.write("%s %s %s %s %s"%(stringsList[i][0], coordSecondMol[i][0], coordSecondMol[i][1], coordSecondMol[i][2], stringsList[i][1])+'\n')
            f.flush()
            os.fsync(f.fileno())
    return fnamearray


def ConvertProbeDimerXYZToTinkerXYZ(poltype,inputxyz,tttxyz,outputxyz,waterbool,probeatoms):
    temp=open(inputxyz,'r')
    resultsinputxyz=temp.readlines()
    temp.close()
    todelete=[]
    for lineidx in range(len(resultsinputxyz)):
        line=resultsinputxyz[lineidx]
        linesplit=line.split()
        if len(linesplit)==0:
            todelete.append(lineidx)
    for index in todelete:
        del resultsinputxyz[index]
 
    temp=open(tttxyz,'r')
    resultstttxyz=temp.readlines()
    temp.close()
    outfile=open(outputxyz,'w')
    outarray=[]
    for lineidx in range(len(resultstttxyz)):
        inputxyzline=resultsinputxyz[lineidx]
        tttxyzline=resultstttxyz[lineidx]
        if lineidx==0:
            totatoms=inputxyzline
        if lineidx!=0:
    
            linesplitinput=inputxyzline.split()
            linesplittttxyz=tttxyzline.split()
            output=''
            output+=' '.join(linesplittttxyz[0:2])+' '
            output+=' '.join(linesplitinput[1:])+' '
            output+=' '.join(linesplittttxyz[5:])
            outarray.append(output)
    
    outfile.write(totatoms)
    for out in outarray:
        outfile.write(out+'\n')
    
    
    totalatomnum=int(totatoms.strip().replace('\n',''))

    if waterbool==True:
        for lineidx in range(len(resultsinputxyz)):
            inputxyzline=resultsinputxyz[lineidx]
            if lineidx==totalatomnum-2:
                outline = str(totalatomnum-2)+ ' ' +inputxyzline.split()[0]+ " " + ' '.join(inputxyzline.split()[1:5]) + " 349 " + str(totalatomnum-1) + " " + str(totalatomnum) + "" 
                outfile.write(outline+'\n')
            elif lineidx==totalatomnum-1:
                outline = str(totalatomnum-1) + ' ' +inputxyzline.split()[0]+ " " + ' '.join(inputxyzline.split()[1:5]) + " 350 " + str(totalatomnum-2)+ ""
                outfile.write(outline+'\n')
            elif lineidx==totalatomnum:
                outline = str(totalatomnum) + ' ' +inputxyzline.split()[0]+ " " + ' '.join(inputxyzline.split()[1:5]) + " 350 " + str(totalatomnum-2)+ ""
                outfile.write(outline+'\n')
    else:
        count=0
        newoutarray=[]
        for lineidx in range(len(resultsinputxyz)):
            inputxyzline=resultsinputxyz[lineidx]
            if lineidx>probeatoms:
                tttindex=lineidx-probeatoms
            else:
                tttindex=lineidx        
            tttxyzline=resultstttxyz[tttindex]
            if lineidx!=0:
                count+=1 
                if count>probeatoms:
                    linesplitinput=inputxyzline.split()
                    linesplittttxyz=tttxyzline.split()
                    output=''
                    atomindex=int(linesplittttxyz[0])
                    linesplittttxyz[0]=str(atomindex+probeatoms)
                    output+=' '.join(linesplittttxyz[0:2])+' '
                    output+=' '.join(linesplitinput[1:])+' '
                    for i in range(6,len(linesplittttxyz)):
                        currentindex=str(int(linesplittttxyz[i])+probeatoms)
                        linesplittttxyz[i]=currentindex
                    output+=' '.join(linesplittttxyz[5:])
                    newoutarray.append(output)
    
        for out in newoutarray:
            outfile.write(out+'\n')

    outfile.close()

def CheckIfLogFileUsingGaussian(poltype,f):
    use_gaus=False
    temp=open(f,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'Entering Gaussian System' in line:
            use_gaus=True
            break 
    return use_gaus



def ReadCounterPoiseAndWriteQMData(poltype,logfilelist):
    temp=open('QM_DATA','w')
    for f in logfilelist:
        use_gaus=False
        use_gaus=CheckIfLogFileUsingGaussian(poltype,f)
        tmpfh=open(f,'r')
        if use_gaus==True:
            for line in tmpfh:
                if 'Counterpoise: corrected energy =' in line or 'Counterpoise corrected energy' in line:
                    dimerenergy=float(line.split()[4])*poltype.Hartree2kcal_mol
            interenergy=dimerenergy
        else:
            for line in tmpfh:
                if 'CP Energy =' in line and 'print' not in line:
                    linesplit=line.split()
                    interenergy=float(linesplit[3])*poltype.Hartree2kcal_mol

        tmpfh.close()
        if interenergy!=None:
            temp.write(f.replace('_sp.log','.xyz')+' '+str(interenergy)+'\n')

    temp.close()




def PlotEnergyVsDistance(poltype,distarray,prefix,radii,depths,reds,vdwtypesarray,count):
    vdwtypes=[str(i) for i in vdwtypesarray]
    vdwtypestring=','.join(vdwtypes)
    qmenergyarray = NormalizeTarget(poltype,'QM_DATA',prefix)
    energyarray = NormalizeTarget(poltype,'SP.dat',prefix)

    energyarray,qmenergyarray,distarray=ScreenHighEnergyPoints(poltype,energyarray,qmenergyarray,distarray)

    def RMSD(c):
        return np.sqrt(np.mean(np.square(np.add(np.subtract(np.array(energyarray),np.array(qmenergyarray)),c))))
    
    r_squared = round(coefficient_of_determination(energyarray,qmenergyarray),2)
    result=fmin(RMSD,.5)
    minRMSD=round(RMSD(result[0]),2)
    if count>1:
        suffix='_boltzman.png'
    else:
        suffix='.png'
    plotname='EnergyVsDistance-'+prefix+'_'+vdwtypestring+suffix
    fig = plt.figure()
    title=prefix+' VdwTypes = '+vdwtypestring
    plt.title(title)
    string=''
    for i in range(len(radii)):
        rad=radii[i]
        depth=depths[i]
        red=reds[i]
        cur='Radius=%s, Depth=%s ,Red=%s'%(round(rad,2),round(depth,2),round(red,2))
        string+=cur
    plt.plot(distarray,energyarray,'bo',label='MM ,'+string)
    plt.plot(distarray,qmenergyarray,'ro',label='QM')
    xpoints=np.array([distarray[i] for i in range(len(distarray))])
    x_new = np.linspace(xpoints.min(),xpoints.max(),500)
    f = interp1d(xpoints,energyarray, kind='quadratic')
    y_smooth=f(x_new)
    plt.plot(x_new,y_smooth,color='blue')
    f = interp1d(xpoints,qmenergyarray, kind='quadratic')
    y_smooth=f(x_new)
    plt.plot(x_new,y_smooth,color='red')
    plt.plot()
    plt.ylabel('Energy (kcal/mol)')
    plt.xlabel('Distance Angstrom '+'RMSD=%s, R^2=%s'%(minRMSD,r_squared))
    plt.legend(loc='best')
    fig.savefig(plotname)


def ReadIntermolecularEnergyMM(poltype,filename):
    energy=None
    with open(filename,'r') as f:
        results=f.readlines()
        for line in results:
            if "Intermolecular Energy :" in line:
                linesplit=line.split()
                energy=linesplit[3]
    return energy

def ReadAnalyzeEnergiesWriteOut(poltype,filenamelist):
    temp=open('SP.dat','w')
    for filename in filenamelist:
        energy=ReadIntermolecularEnergyMM(poltype,filename)
        temp.write(filename+' '+str(energy)+'\n')
    temp.close() 


def WriteInitialPrmFile(poltype,vdwtypesarray,initialradii,initialdepths,minradii,maxradii,mindepths,maxdepths,initialreds,minreds,maxreds):

    temp=open('INITIAL.PRM','w')
    for i in range(len(vdwtypesarray)):
        vdwtype=str(vdwtypesarray[i])
        initialradius=str(initialradii[i])
        initialdepth=str(initialdepths[i])
        initialred=str(initialreds[i])
        minradius=str(minradii[i])
        maxradius=str(maxradii[i])
        mindepth=str(mindepths[i])
        maxdepth=str(maxdepths[i])
        minred=str(minreds[i])
        maxred=str(maxreds[i])
        line='vdw '+vdwtype+' '+initialradius+' '+initialdepth+' '+initialred+' '+minradius+' '+maxradius+' '+mindepth+' '+maxdepth+' '+minred+' '+maxred+'\n'
        temp.write(line)

    temp.close()

def readTXYZ(poltype,TXYZ):
    temp=open(TXYZ,'r')
    lines = temp.readlines()[1:] #TINKER coordinate starts from second line
    atoms=[];coord=[]
    order=[];types=[];connections=[]
    for line in lines:
        data=line.split()
        order.append(data[0])
        types.append(data[5])
        connections.append(data[6:])
        atoms.append(data[1])
        coord.append([float(data[2]), float(data[3]), float(data[4])])
    return atoms,coord,order, types, connections

def TXYZ2COM(poltype,TXYZ,comfname,chkname,maxdisk,maxmem,numproc,mol,probeatoms):
    data = readTXYZ(poltype,TXYZ)
    atoms = data[0];coord = data[1]
    opt.write_com_header(poltype,comfname,chkname,maxdisk,maxmem,numproc)
    tmpfh = open(comfname, "a")
    if ('I ' in poltype.mol.GetSpacedFormula()):
        poltype.espbasisset='gen'
        iodinebasissetfile=poltype.iodineespbasissetfile 
        basissetfile=poltype.espbasissetfile 
    
    opstr="#P %s/%s Sp Counterpoise=2" % (poltype.espmethod,poltype.espbasisset)

    if ('I ' in poltype.mol.GetSpacedFormula()):
        opstr+=' pseudo=read'
    string=' MaxDisk=%s \n'%(maxdisk)
    opstr+=string
    mul=mol.GetTotalSpinMultiplicity()
    chg=mol.GetTotalCharge()
    tmpfh.write(opstr)
    commentstr = poltype.molecprefix + " Gaussian SP Calculation on " + gethostname()
    tmpfh.write('\n%s\n\n' % commentstr)
    tmpfh.write('%d %d %d %d %d %d\n' % (chg,mul,chg,mul,0,1))
    fragelements=[]
    for n in range(len(atoms)):
        if n>=len(atoms)-probeatoms:
            ele=atoms[n]
            if ele not in fragelements:
                 fragelements.append(ele)

            tmpfh.write("%3s%s             %14.7f%14.7f%14.7f\n"%(atoms[n],'(Fragment=2)',float(coord[n][0]),float(coord[n][1]),float(coord[n][2]))) 
        else:
            tmpfh.write("%3s%s             %14.7f%14.7f%14.7f\n"%(atoms[n],'(Fragment=1)',float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))   
    tmpfh.write("\n")

    if ('I ' in poltype.mol.GetSpacedFormula()):
        formulalist=poltype.mol.GetSpacedFormula().lstrip().rstrip().split()
        elementtobasissetlines=GenerateElementToBasisSetLines(poltype,poltype.basissetpath+basissetfile)
        for element,basissetlines in elementtobasissetlines.items():
            if element in poltype.mol.GetSpacedFormula() or element in fragelements:
                for line in basissetlines: 
                    tmpfh.write(line)


        temp=open(poltype.basissetpath+iodinebasissetfile,'r')
        results=temp.readlines()
        temp.close()
        for line in results:
            if '!' not in line:
                tmpfh.write(line)
        tmpfh.write('\n')

        

    tmpfh.write("\n")
    tmpfh.close()

def GenerateElementToBasisSetLines(poltype,basissetfile):
    elementtobasissetlines={}
    temp=open(basissetfile,'r')
    results=temp.readlines()
    temp.close()
    lines=[]
    element=None
    for line in results:
        linesplit=line.split()
        if len(linesplit)==2 and linesplit[0].isalpha() and linesplit[1]=='0':
            if element!=None:
                elementtobasissetlines[element]=lines
            element=linesplit[0]
            lines=[line]
        else:
            lines.append(line)
    return elementtobasissetlines
 

def ReadInBasisSet(poltype,tmpfh,normalelementbasissetfile,otherelementbasissetfile,space):
    newtemp=open(poltype.basissetpath+normalelementbasissetfile,'r')
    results=newtemp.readlines()
    newtemp.close()
    for line in results:
        if '!' not in line:
            tmpfh.write(space+line)


    newtemp=open(poltype.basissetpath+otherelementbasissetfile,'r')
    results=newtemp.readlines()
    newtemp.close()
    for line in results:
        if '!' not in line:
            tmpfh.write(space+line)
    return tmpfh




def CreatePsi4SPInputFile(poltype,TXYZ,mol,maxdisk,maxmem,numproc,probeatoms):
    data = readTXYZ(poltype,TXYZ)
    atoms = data[0];coord = data[1]
    mul=mol.GetTotalSpinMultiplicity()
    chg=mol.GetTotalCharge()
    inputname=TXYZ.replace('.xyz','_sp.psi4')

    temp=open(inputname,'w')
    temp.write('molecule { '+'\n')
    temp.write('%d %d\n' % (chg, mul))
    for n in range(len(atoms)):
        if n==len(atoms)-probeatoms:
            temp.write('--'+'\n')
            temp.write('%d %d\n' % (0, 1))
        if n>=len(atoms)-probeatoms:
            temp.write("%3s             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2]))) 
        else:
            temp.write("%3s             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))    

    temp.write('}'+'\n')
    temp.write('memory '+maxmem+'\n')
    temp.write('set_num_threads(%s)'%(numproc)+'\n')
    temp.write('psi4_io.set_default_path("%s")'%(poltype.scrtmpdirpsi4)+'\n')
    temp.write('set maxiter '+str(poltype.scfmaxiter)+'\n')
    temp.write('set freeze_core True'+'\n')
    temp.write('set PROPERTIES_ORIGIN ["COM"]'+'\n')
    if poltype.allowradicals==True:
        temp.write('set reference uhf '+'\n')

    spacedformulastr=mol.GetSpacedFormula()
    if ('I ' in spacedformulastr):
        temp.write('basis {'+'\n')
        temp.write('['+' '+poltype.espbasissetfile+' '+poltype.iodineespbasissetfile +' '+ ']'+'\n')
        temp=ReadInBasisSet(poltype,temp,poltype.espbasissetfile,poltype.iodineespbasissetfile,'')
        temp.write('}'+'\n')
        temp.write("e_dim= energy('%s',bsse_type='cp')" % (poltype.espmethod.lower())+'\n')
    else:
        temp.write("e_dim= energy('%s/%s',bsse_type='cp')" % (poltype.espmethod.lower(),poltype.espbasisset)+'\n')
    temp.write('\n')
    temp.write("psi4.print_out('CP Energy = %10.6f' % (e_dim))"+'\n')
    temp.write('clean()'+'\n')
    temp.close()
    temp=open(inputname,'r')
    results=temp.readlines()
    temp.close()
    outputname=os.path.splitext(inputname)[0] + '.log'
    return inputname,outputname



def WriteOutCartesianXYZ(poltype,mol,filename):
    output=open(filename,'w')
    atomcounter=0
    coordarray=[]
    atomiter=openbabel.OBMolAtomIter(mol)
    an = pyasl.AtomicNo()
    for atom in atomiter:
        atomcounter+=1
        atomsymb=an.getElSymbol(atom.GetAtomicNum())
        x=str(atom.GetX())
        y=str(atom.GetY())
        z=str(atom.GetZ())
        xyzline=atomsymb+' '+x+' '+y+' '+z
        coordarray.append(xyzline)
            
    output.write(str(atomcounter)+'\n')
    for coord in coordarray:
        output.write(coord+'\n')
    output.close()


def GenerateSPInputFiles(poltype,filenamearray,mol,probeatoms):
    qmfilenamearray=[]
    for filename in filenamearray:
        if poltype.use_gaus==True:
            qmfilename=filename.replace('.xyz','_sp.com')
            chkfilename=filename.replace('.xyz','_sp.chk')
            TXYZ2COM(poltype,filename,qmfilename,chkfilename,poltype.maxdisk,poltype.maxmem,poltype.numproc,mol,probeatoms)
        else:
            qmfilename,outputname=CreatePsi4SPInputFile(poltype,filename,mol,poltype.maxdisk,poltype.maxmem,poltype.numproc,probeatoms)
        qmfilenamearray.append(qmfilename)
    return qmfilenamearray


def ExecuteSPJobs(poltype,qmfilenamearray,prefix):
    jobtooutputlog={}
    listofjobs=[]
    fulljobtooutputlog={}
    jobtoinputfilepaths={}
    jobtooutputfiles={}
    jobtoabsolutebinpath={}
    outputfilenames=[]
    for i in range(len(qmfilenamearray)):
        filename=qmfilenamearray[i]
        if poltype.use_gaus==True:
            cmdstr = 'GAUSS_SCRDIR='+poltype.scrtmpdirgau+' '+poltype.gausexe+' '+filename
            outputname=filename.replace('.com','.log')
            executable=poltype.gausexe
            scratchdir=poltype.scrtmpdirgau
        else:
            outputname=filename.replace('.psi4','.log')
            cmdstr='psi4 '+filename+' '+outputname
            executable='psi4'
            scratchdir=poltype.scrtmpdirpsi4
        abspath=poltype.which(executable)
        inputfilepath=os.path.join(os.getcwd(),filename) 
        finished,error=poltype.CheckNormalTermination(outputname)
        if finished==True and error==False:
            pass
        else:
            listofjobs.append(cmdstr)
            jobtooutputlog[cmdstr]=os.getcwd()+r'/'+outputname
            jobtoinputfilepaths[cmdstr]=[inputfilepath]
            jobtooutputfiles[cmdstr]=[outputname]
            jobtoabsolutebinpath[cmdstr]=abspath

        outputfilenames.append(outputname)
    lognames=[]
    for job in listofjobs:
        log=jobtooutputlog[job]
        lognames.append(os.path.abspath(poltype.logfname))
    fulljobtolog=dict(zip(listofjobs, lognames)) 
    fulljobtooutputlog.update(jobtooutputlog)
    jobtologlistfilenameprefix=os.getcwd()+r'/'+'QMSPJobToLog'+'_'+prefix
    if poltype.externalapi!=None:
        if len(listofjobs)!=0:
            call.CallExternalAPI(poltype,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,scratchdir,jobtologlistfilenameprefix)
        finshedjobs,errorjobs=poltype.WaitForTermination(fulljobtooutputlog,False)
    else:
        finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(fulljobtooutputlog,True)

    return outputfilenames


def GrabVdwParameters(poltype,vdwtype):
    CheckAllVdwTypesExist(poltype,poltype.key4fname)
    temp=open(poltype.key4fname,'r')
    results=temp.readlines()
    temp.close()
    if 'S' in poltype.vdwprmtypestofit:
        radbound=.1
    else:
        radbound=0
    if 'T' in poltype.vdwprmtypestofit:
        depthbound=.1
    else:
        depthbound=0
    for line in results:
        if 'vdw' in line:
            if str(vdwtype) in line:
                if str(vdwtype) in poltype.fixvdwtyperadii:
                    rbound=0
                else:
                    rbound=radbound
                linesplit=line.split()
                radius=float(linesplit[2])
                minvdwradius=radius-rbound*radius
                maxvdwradius=radius+rbound*radius
                depth=float(linesplit[3])
                minvdwdepth=depth-depthbound*depth
                maxvdwdepth=depth+depthbound*depth
                if depth>maxvdwdepth:
                    depth=maxvdwdepth-.001
                if len(linesplit)==5 and poltype.fitred==True:
                    red=float(linesplit[4])
                else:
                    red=1
                minred=red-.1*red
                maxred=1
                return radius,depth,minvdwradius,maxvdwradius,minvdwdepth,maxvdwdepth,red,minred,maxred


def GenerateReferenceDistances(indextoreferencecoordinate,indextomolecule,indextotargetatom,indextoreferenceelement,vdwradius):
    indexpairtoreferencedistance={}
    indexpairtobounds={}
    indices=list(indextoreferencecoordinate.keys())
    allpairs=list(itertools.combinations(indices, 2)) 
    tol=1
    distpairs=[]
    for pair in allpairs:
        index1=pair[0]
        index2=pair[1]
        molecule1=indextomolecule[index1]
        molecule2=indextomolecule[index2]
        coord1=indextoreferencecoordinate[index1]
        coord2=indextoreferencecoordinate[index2]
        target1=indextotargetatom[index1]
        target2=indextotargetatom[index2]
        element1=indextoreferenceelement[index1]
        element2=indextoreferenceelement[index2]
        vdwradius1=vdwradius[element1]
        vdwradius2=vdwradius[element2]
        if molecule1==molecule2:
            targetdistance=np.linalg.norm(coord1-coord2)
            bound=[targetdistance,targetdistance]
        else:
            if target1=='target' and target2=='target': # then use vdw distances
                targetdistance=(vdwradius1**3+vdwradius2**3)/(vdwradius1**2+vdwradius2**2)
                if targetdistance<2:
                    targetdistance=2
                bound=[targetdistance,targetdistance]
                distpairs.append(pair)
        indexpairtoreferencedistance[tuple(pair)]=targetdistance 
        indexpairtobounds[tuple(pair)]=bound

    for pair in allpairs:
        index1=pair[0]
        index2=pair[1]
        molecule1=indextomolecule[index1]
        molecule2=indextomolecule[index2]
        coord1=indextoreferencecoordinate[index1]
        coord2=indextoreferencecoordinate[index2]
        target1=indextotargetatom[index1]
        target2=indextotargetatom[index2]
        element1=indextoreferenceelement[index1]
        element2=indextoreferenceelement[index2]
        vdwradius1=vdwradius[element1]
        vdwradius2=vdwradius[element2]
        if molecule1==molecule2:
            targetdistance=np.linalg.norm(coord1-coord2)
            bound=[targetdistance,targetdistance]
        else:
            if target1=='target' and target2=='target': # then use vdw distances
                targetdistance=(vdwradius1**3+vdwradius2**3)/(vdwradius1**2+vdwradius2**2)
                if targetdistance<2:
                    targetdistance=2

                bound=[targetdistance-tol,targetdistance+tol]
            else:
                if targetdistance<2:
                    targerdistance=2
                bound=[targetdistance,100] # high upper bound, large range to not add cost in cost function                

        indexpairtoreferencedistance[tuple(pair)]=targetdistance 
        indexpairtobounds[tuple(pair)]=bound

 
    return indexpairtoreferencedistance,indexpairtobounds,distpairs


def GrabPairwiseDistance(pair,indexpairtoreferencedistance):
    if pair in indexpairtoreferencedistance.keys():
        distance=indexpairtoreferencedistance[pair]
    elif pair[::-1] in indexpairtoreferencedistance.keys():
        distance=indexpairtoreferencedistance[pair[::-1]]
    return distance



def GenerateInitialDictionaries(coords1,coords2,atoms1,atoms2,p1,p2):
    indextoreferencecoordinate={}
    indextoreferenceelement={}
    indextomolecule={}
    indextotargetatom={}
    count=0
    for i in range(len(coords1)):
        coords=np.array(coords1[i])
        element=atoms1[i]

        indextoreferencecoordinate[count]=coords
        indextoreferenceelement[count]=element
        indextomolecule[count]='molecule1'
        if count==p1:
            string='target'
        else:
            string='steric'
        indextotargetatom[count]=string
        count+=1     
    p2=p2+len(coords1)
    for i in range(len(coords2)):
        coords=np.array(coords2[i])
        element=atoms2[i]
        indextoreferencecoordinate[count]=coords
        indextoreferenceelement[count]=element
        indextomolecule[count]='molecule2'
        if count==p2:
            string='target'
        else:
            string='steric'
        indextotargetatom[count]=string

        count+=1     
    return indextoreferencecoordinate,indextoreferenceelement,indextomolecule,indextotargetatom


def Generate3DMeshGrid(poltype,moleculecoords,probecoords):
    grid=[]
    meshsize=.5 # angstroms, half smallest bond length
    maxdistmolecule,maxcoordmolecule=FindLongestDistanceInMolecule(poltype,moleculecoords)
    maxdistprobe,maxcoordprobe=FindLongestDistanceInMolecule(poltype,probecoords)
    sidebuffer=2*(maxdistprobe+2) # number of angstroms away from end of molecule on both sides, so if water molecule is probe, longest length of water +2 ang, then x 2, add this to both sides of longest length of molecule
    for i in range(len(maxcoordmolecule)):
        coord=maxcoordmolecule[i]
        lower=coord-sidebuffer
        upper=coord+maxdistmolecule+sidebuffer
        row=np.arange(lower,upper,meshsize)
        grid.append(row)
    xv, yv , zv= np.meshgrid(grid[0], grid[1],grid[2])
    return xv,yv,zv


def FindLongestDistanceInMolecule(poltype,coords):
    veclist=[]
    for i in range(len(coords)):
        vec=coords[i]
        veclist.append(vec)
    pairs=list(itertools.combinations(veclist, 2))
    disttopair={}
    for pairidx in range(len(pairs)):
        pair=pairs[pairidx]
        dist=np.linalg.norm(np.array(pair[0])-np.array(pair[1]))
        disttopair[dist]=pair
    distances=list(disttopair.keys())
    maxdist=max(distances)
    maxpair=disttopair[maxdist]
    maxcoord=maxpair[0] # just choose first one as starting point
    return maxdist,maxcoord

def GrabMoleculeCoords(poltype,indextoreferencecoordinate,indextomolecule,moleculekey):
    moleculeatoms=[]
    for i in indextoreferencecoordinate.keys():
        molecule=indextomolecule[i]
        if molecule==moleculekey:
            moleculeatoms.append(i)
    coords=[indextoreferencecoordinate[i] for i in moleculeatoms]
    return coords,moleculeatoms

def FindBestGridPoints(poltype,x,y,z,moleculecoords,probecoords,indextoreferencecoordinate,refcoord,refdistance):
    bestgridpoints=[]
    maxdistprobe,maxcoordprobe=FindLongestDistanceInMolecule(poltype,probecoords)
    differencetopoint={}
    for i in range(len(x)):
        for j in range(len(x[0])):
            for k in range(len(x[0,0])):
                xvalue=x[i,j,k]
                yvalue=y[i,j,k]
                zvalue=z[i,j,k]
                point=np.array([xvalue,yvalue,zvalue])
                nosteric=True
                for coord in moleculecoords:
                    if not np.array_equal(np.array(refcoord),np.array(coord)):
                        distance=np.linalg.norm(np.array(coord)-np.array(point))
                        if distance<refdistance: # want to be at least refdistance away from other atoms and also target atom, if any closer then not probing the target atom  
                           nosteric=False
                           break
                if nosteric==True:
                    distance=np.linalg.norm(np.array(refcoord)-np.array(point))-refdistance
                    difference=np.abs(distance) 
                    differencetopoint[difference]=point

    sorted_dict = {}
    sorted_keys = sorted(differencetopoint.keys())  
    for w in sorted_keys:
        sorted_dict[w] = differencetopoint[w]
    count=0
    diffs=[]
    for key,value in sorted_dict.items():
        count+=1
        bestgridpoints.append(value)
        diffs.append(key)
        if count>poltype.vdwmaxtinkergridpoints:
            break
    return bestgridpoints




def TranslateProbe(poltype,indextoreferencecoordinate,bestgridpoint,probecoords,probeatoms,probedonorindex):
    transvec=np.array(bestgridpoint)-np.array(indextoreferencecoordinate[probedonorindex]) 
    for i in range(len(probecoords)): 
        probecoord=probecoords[i]
        probeindex=probeatoms[i]
        newcoord=np.array(probecoord)+transvec
        indextoreferencecoordinate[probeindex]=newcoord

    return indextoreferencecoordinate


def GenerateStartingDimers(poltype,indextoreferencecoordinate,indexpairtoreferencedistance,distpairs,indextomolecule):
    moleculecoords,moleculeatoms=GrabMoleculeCoords(poltype,indextoreferencecoordinate,indextomolecule,'molecule1')
    probecoords,probeatoms=GrabMoleculeCoords(poltype,indextoreferencecoordinate,indextomolecule,'molecule2')
    x,y,z=Generate3DMeshGrid(poltype,moleculecoords,probecoords) 
    distpair=distpairs[0]
    for i in range(len(distpair)):
        index=distpair[i]
        molecule=indextomolecule[index]
        if molecule=='molecule1':
            refcoord=indextoreferencecoordinate[index]
            refdistance=indexpairtoreferencedistance[distpair]
        elif molecule=='molecule2':
            probeindex=index
    bestgridpoints=FindBestGridPoints(poltype,x,y,z,moleculecoords,probecoords,indextoreferencecoordinate,refcoord,refdistance)
    listofindextoreferencecoordinate=[]
    for bestgridpoint in bestgridpoints:
        newindextoreferencecoordinate=indextoreferencecoordinate.copy()
        newindextoreferencecoordinate=TranslateProbe(poltype,newindextoreferencecoordinate,bestgridpoint,probecoords,probeatoms,probeindex)
        listofindextoreferencecoordinate.append(newindextoreferencecoordinate)


    return listofindextoreferencecoordinate





def optimizedimers(poltype,atoms1, atoms2, coords1, coords2, p1, p2, oridimer,vdwradius,mol,probemol,probeidxtosymclass):
    indextoreferencecoordinate,indextoreferenceelement,indextomolecule,indextotargetatom=GenerateInitialDictionaries(coords1,coords2,atoms1,atoms2,p1,p2) 
    indexpairtoreferencedistance,indexpairtobounds,distpairs=GenerateReferenceDistances(indextoreferencecoordinate,indextomolecule,indextotargetatom,indextoreferenceelement,vdwradius)
    indexpairtoreferencedistanceoriginal=indexpairtoreferencedistance.copy()
    listofindextoreferencecoordinate=GenerateStartingDimers(poltype,indextoreferencecoordinate,indexpairtoreferencedistance,distpairs,indextomolecule)
    minstructs=[]
    for i in range(len(listofindextoreferencecoordinate)):
        indextoreferencecoordinate=listofindextoreferencecoordinate[i]
        dimer=oridimer.replace('.xyz','_'+str(i)+'.xyz')
        with open(dimer, "w") as f:
          f.write(str(len(indextoreferencecoordinate.keys()))+"\n")
          f.write("\n")
          for index,coordinate in indextoreferencecoordinate.items():
              element=indextoreferenceelement[index]
              f.write("%3s %12.5f%12.5f%12.5f\n"%(element, coordinate[0], coordinate[1], coordinate[2]))
        if 'water' in dimer:
            waterbool=True
        else:
            waterbool=False
        outputxyz=dimer.replace('.xyz','-tinkermin.xyz')
        probeatoms=len(atoms2)
        ConvertProbeDimerXYZToTinkerXYZ(poltype,dimer,poltype.xyzoutfile,outputxyz,waterbool,probeatoms)
        outputxyz=MinimizeDimer(poltype,outputxyz,poltype.key4fname)
        if outputxyz!=None:
            minstructs.append(outputxyz)
    return minstructs

def MinimizeDimer(poltype,inputxyz,keyfile):
    torminlogfname=inputxyz.replace('.xyz','.out')
    mincmdstr = poltype.minimizeexe+' '+inputxyz+' -k '+keyfile+' 0.1'+' '+'>'+torminlogfname
    term,error=poltype.CheckNormalTermination(torminlogfname)
    errorjobs=[]
    if term==True and error==False:
        pass
    else:
        poltype.call_subsystem([mincmdstr],False)
        temp={mincmdstr:torminlogfname} 
        finishedjobs,errorjobs=poltype.WaitForTermination(temp,True)
    if len(errorjobs)==0:
        finaloutputxyz=inputxyz+'_2'
    else:
        finaloutputxyz=None
    return finaloutputxyz


def ConvertTinktoXYZ(poltype,filename,newfilename):
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
            newline=linesplit[1]+' '+linesplit[2]+' '+linesplit[3]+' '+linesplit[4]+'\n'
            tempwrite.write(newline)
            tempwrite.flush()
            os.fsync(tempwrite.fileno())
    temp.close()
    tempwrite.close()





def GenerateProbePathNames(poltype,vdwprobenames,moleculexyz):
    probepaths=[]
    if 'water' in vdwprobenames:
        path=poltype.vdwprobepathname+'water.xyz'
        probepaths.append(path)
    if poltype.homodimers==True:
        probepaths.append(moleculexyz)
    return probepaths


def readXYZ(poltype,xyz):
  atoms  = np.loadtxt(xyz, usecols=(0,), dtype='str', unpack=True, skiprows=2)
  coords = np.loadtxt(xyz, usecols=(1,2,3), dtype='float', unpack=False, skiprows=2)
  return atoms,coords


def CheckBuriedAtoms(poltype,indexarray,molecule,zeroindex=False):
    indicestodelete=[]
    for i in range(len(indexarray)):
        index=indexarray[i]
        if zeroindex==False:
            idx=index
        else:
            idx=index+1
        atom=molecule.GetAtom(idx)
        valence=atom.GetExplicitValence()
        hyb=atom.GetHyb() 
        if valence==4 and hyb==3:
            indicestodelete.append(i)
    for index in indicestodelete:
        del indexarray[index]

    return indexarray


def GrabWaterTypes(poltype):
    probeidxtosymclass={1:349,2:350,3:350}
    return probeidxtosymclass


def GenerateInitialProbeStructures(poltype,missingvdwatomindices):
    molecules = [poltype.xyzfname]
    probes = GenerateProbePathNames(poltype,poltype.vdwprobenames,poltype.xyzfname)
    vdwradius = {"H" : 1.20, "Li": 1.82, "Na": 2.27, "K": 2.75, "Rb": 3.03, "Cs": 3.43, \
                 "Be": 1.53, "Mg": 1.73, "Ca": 2.31, "B": 1.92, "C": 1.70, "N": 1.55, "O":1.52, \
                 "P" : 1.80, "S" : 1.80, "F" : 1.47, "Cl":1.75, "Br":1.85, "Zn":1.39,'I':4.61}           
    dimernames=[]
    probeindices=[]
    moleculeindices=[]
    numberprobeatoms=[]
    for mol in molecules:
        atoms1,coords1,order1, types1, connections1=readTXYZ(poltype,mol)
        mol_spots = missingvdwatomindices.copy()
        mol_spots = [i-1 for i in mol_spots]
        moldimernames=[]
        for prob in probes:
            probename=os.path.basename(prob)
            prefix=poltype.molstructfname.split('.')[0]
            if prefix in prob:
                probemol=poltype.mol
                atoms2=atoms1.copy()
                coords2=coords1.copy()
            else:
                atoms2,coords2=readXYZ(poltype,prob)
                probemol=opt.load_structfile(poltype,prob)
            probeidxtosymclass,symmetryclass=symm.gen_canonicallabels(poltype,probemol) 
            prob_spots=[]

            for symclass in list(probeidxtosymclass.values()): 
                keys= GrabKeysFromValue(poltype,probeidxtosymclass,symclass)       
                probeidx=keys[0]-1 # shift to 0 index
                if probeidx not in prob_spots:
                    prob_spots.append(probeidx)
            addprobetype=False
            if 'water' in prob:
                probeidxtosymclass=GrabWaterTypes(poltype)
            else:
                addprobetype=True
            prob_spots = CheckBuriedAtoms(poltype,prob_spots,probemol,zeroindex=True)
            alreadyused=[]
            for p1 in mol_spots:
                probelist=[]
                probeindiceslist=[]
                moleculeindiceslist=[]
                moltype=poltype.idxtosymclass[p1+1]
                newls=[moltype]
                for p2 in prob_spots:
                    if addprobetype==True:
                        probetype=probeidxtosymclass[p2+1]
                        newls.append(probetype)
                    setls=set(newls)
                    if setls in alreadyused and len(setls)>1:
                        continue
                    else:
                        alreadyused.append(setls)
                        e1=atoms1[p1]
                        e2=atoms2[p2]
                        if (e1=='H' and e2=='H') or (e1=='O' and e2=='O'):
                            continue
                        dimer = mol[:-4] + "-" + probename[:-4] + "_" + str("%d_%d"%(p1+1,len(atoms1)+p2+1)) + ".xyz"
                        minstructs=optimizedimers(poltype,atoms1, atoms2, coords1, coords2, p1, p2, dimer,vdwradius,poltype.mol,probemol,probeidxtosymclass)
                        probelist.append(minstructs)
                        probevalue=p2+1+len(atoms1)
                        ls=[]
                        for k in range(len(minstructs)):
                            ls.append(probevalue)
                        probeindiceslist.append(ls)
                        molvalue=p1+1
                        ls=[]
                        for k in range(len(minstructs)):
                            ls.append(molvalue)
                        moleculeindiceslist.append(ls)
                moleculeindices.append(moleculeindiceslist)
                probeindices.append(probeindiceslist)
                moldimernames.append(probelist)
                numberprobeatoms.append(probemol.NumAtoms())
    return moldimernames,probeindices,moleculeindices,numberprobeatoms




def GrabKeysFromValue(poltype,dic,thevalue):
    keylist=[]
    for key,value in dic.items():
        if value==thevalue:
            keylist.append(key)
    return keylist


def ReplaceParameterFileHeader(poltype,paramhead,keyfile):
    tempname='temp.key'
    temp=open(keyfile,'r')
    results=temp.readlines()
    temp.close() 
    temp=open(tempname,'w')
    for line in results:
        if 'parameter' in line and '#' not in line:
            newline='parameters '+paramhead+'\n'
        else:
            newline=line
        temp.write(newline)
    temp.close() 
    os.remove(keyfile)
    os.rename(tempname,keyfile)
  
def CheckIfFittingCompleted(poltype,prefix):
    check=False
    files=os.listdir()
    for f in files:
        if prefix in f and '.png' in f:
            check=True
            break
    return check


def intersection(poltype,lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def CombineProbesThatNeedToBeFitTogether(poltype,probeindices,moleculeindices,fullprefixarrays,fulldistarrays,alloutputfilenames):
    newprobeindices=[]
    newmoleculeindices=[]
    newprefixarrays=[]
    newdistarrays=[]
    newoutputfilenames=[]
    indextoindexlist={}
    molindextotypenum={}
    group=[]
    for i in range(len(moleculeindices)):
        moleculelist=moleculeindices[i]
        index=moleculelist[0]
        typenum=poltype.idxtosymclass[index]
        checkboth=False
        filename=fullprefixarrays[i][0]
        if 'water' not in filename:
            checkboth=True
        probelist=probeindices[i]
        probeindex=probelist[0]
        molindextotypenum[index]=typenum
        ls=[index]
        if checkboth==True:
            newindex=int(probeindex-len(poltype.idxtosymclass.keys()))
            typenum=poltype.idxtosymclass[newindex]   
            molindextotypenum[probeindex]=typenum
            ls.append(probeindex)
        group.append(ls)
    groupsym=[]
    for i in range(len(group)):
        ls=group[i]
        symls=[molindextotypenum[k] for k in ls]
        groupsym.append(symls)
    indicesgrouped=[]
    for i in range(len(groupsym)):
        symls=groupsym[i]
        ls=group[i]
        for j in range(len(groupsym)):
            if j!=i:
                othersymls=groupsym[j]
                otherls=group[j]
                
                inter=intersection(poltype,symls, othersymls)  
                index=ls[0]
                otherindex=otherls[0]
                if len(inter)!=0:
                    if i not in indextoindexlist.keys():
                        indextoindexlist[i]=[]
                    if i not in indextoindexlist[i]:
                        indextoindexlist[i].append(i)
                    if j not in indextoindexlist[i]:
                        indextoindexlist[i].append(j)
                    if index not in indicesgrouped:
                        indicesgrouped.append(index)
                    if otherindex not in indicesgrouped:
                        indicesgrouped.append(otherindex)
    for i in range(len(group)):
        ls=group[i]
        index=ls[0]
        if index not in indicesgrouped:
            indicesgrouped.append(index)
            if i not in indextoindexlist.keys():
                indextoindexlist[i]=[]
            indextoindexlist[i].append(i)


    indicesalreadydone=[]
    for index,indexlist in indextoindexlist.items():
        if index not in indicesalreadydone:
            tempmoleculeindices=[]
            tempprobeindices=[]
            tempprefixes=[]
            tempdistarrays=[]
            tempfilenames=[]
            for idx in indexlist:
                indicesalreadydone.append(idx)  
                moleculelist=moleculeindices[idx]
                probelist=probeindices[idx]
                prefixes=fullprefixarrays[idx]
                distarrays=fulldistarrays[idx]
                filenames=alloutputfilenames[idx]
                tempmoleculeindices.append(moleculelist)
                tempprobeindices.append(probelist)
                tempprefixes.append(prefixes)
                tempdistarrays.append(distarrays)     
                tempfilenames.append(filenames)
            
            newprobeindices.append(tempprobeindices)
            newmoleculeindices.append(tempmoleculeindices)
            newprefixarrays.append(tempprefixes)
            newdistarrays.append(tempdistarrays)
            newoutputfilenames.append(tempfilenames)  

    return newprobeindices,newmoleculeindices,newprefixarrays,newdistarrays,newoutputfilenames
 

def RemoveIgnoredIndices(poltype,probeindices,moleculeindices,moleculeprobeindicestoignore):
    finalprobeindices=[]
    finalmoleculeindices=[]
    for i in range(len(probeindices)):
        probeindexlist=probeindices[i]
        moleculeindexlist=moleculeindices[i]
        newprobeindexlist=[]
        newmoleculeindexlist=[]
        for j in range(len(probeindexlist)):
            probeindex=probeindexlist[j]
            moleculeindex=moleculeindexlist[j]
            ls=[moleculeindex,probeindex]
            if ls in moleculeprobeindicestoignore:
                pass
            else:
                newmoleculeindexlist.append(moleculeindex)
                newprobeindexlist.append(probeindex)
        if len(newprobeindexlist)!=0:
            finalprobeindices.append(newprobeindexlist)
        if len(newmoleculeindexlist)!=0:
            finalmoleculeindices.append(newmoleculeindexlist)


    return finalprobeindices,finalmoleculeindices


def CheckIfProbeIsTooFarOrTooClose(poltype,mol,probeindex,moleculeindex,probeatoms):
    checktoofar=False
    moleculeatom=mol.GetAtom(moleculeindex)
    probeatom=mol.GetAtom(probeindex)
    moleculeatomcoords=np.array([moleculeatom.GetX(),moleculeatom.GetY(),moleculeatom.GetZ()])
    probeatomcoords=np.array([probeatom.GetX(),probeatom.GetY(),probeatom.GetZ()])
    dist=np.linalg.norm(probeatomcoords-moleculeatomcoords)
    if dist>=5 or dist<1.6:
        checktoofar=True
    atomiter=openbabel.OBMolAtomIter(mol)
    total=0
    for atom in atomiter:
        total+=1
    tol=.5
    maxatomindex=total-probeatoms-1
    atomiter=openbabel.OBMolAtomIter(mol)
    for atom in atomiter:
        atomindex=atom.GetIndex()
        theatomindex=atomindex+1
        if theatomindex!=(moleculeindex) and atomindex<=maxatomindex:
            othermoleculeatomcoords=np.array([atom.GetX(),atom.GetY(),atom.GetZ()])
            otherdist=np.linalg.norm(probeatomcoords-othermoleculeatomcoords)
            if otherdist+tol<dist:
                checktoofar=True
        



    return checktoofar


def LonePairAboveBelowPlane(moleculeatom,mol,moleculeatomcoords,moleculeindex):
    atomiter=openbabel.OBAtomAtomIter(moleculeatom)
    natoms=[]
    found=False
    for natom in atomiter:
        natomisinring=natom.IsInRing()
        if natomisinring==True:
            atomatomiter=openbabel.OBAtomAtomIter(natom)
            for nnatom in atomatomiter:
                nnatomidx=nnatom.GetIdx()
                if nnatomidx!=moleculeindex and nnatom.IsInRing()==True:
                    natoms.append(natom)
                    natoms.append(nnatom)
                    found=True 
                    break
            if found==True:
                break
    npos=[np.array([natom.GetX(),natom.GetY(),natom.GetZ()]) for natom in natoms]
    displacements=[pos-moleculeatomcoords for pos in npos]
    dists=[np.linalg.norm(natomcoords-moleculeatomcoords) for natomcoords in npos]
    lonepairvec=np.zeros(3)
    unitvectors=[]
    for i in range(len(displacements)):
        disp=displacements[i]
        dist=dists[i]
        unitvec=disp/dist
        unitvectors.append(unitvec)
    crossprod=np.cross(unitvectors[0],unitvectors[1])
    lonepairvecs=[crossprod,-1*crossprod]
    return lonepairvecs

def LonePairBisector(moleculeatom,mol,moleculeatomcoords,moleculeindex):
    atomiter=openbabel.OBAtomAtomIter(moleculeatom)
    natoms=[]
    for natom in atomiter:
        natoms.append(natom)
    npos=[np.array([natom.GetX(),natom.GetY(),natom.GetZ()]) for natom in natoms]
    displacements=[pos-moleculeatomcoords for pos in npos]
    dists=[np.linalg.norm(natomcoords-moleculeatomcoords) for natomcoords in npos]
    lonepairvec=np.zeros(3)
    unitvectors=[]
    for i in range(len(displacements)):
        disp=displacements[i]
        dist=dists[i]
        unitvec=disp/dist
        unitvectors.append(unitvec)
    disp=unitvectors[0]-unitvectors[1]
    disp=.5*disp
    lonepairvec=disp-moleculeatomcoords
    norm=np.linalg.norm(lonepairvec)
    lonepairvec=-1*lonepairvec/norm 
    lonepairvecs=[lonepairvec]
    return lonepairvecs


def AddLonePairPoints(poltype,newdimerfiles,newprobeindices,newmoleculeindices,newnumberprobeatoms,moleculeindex,probeindex,mol,probeatoms,newfilename):
    added=False
    moleculeatom=mol.GetAtom(moleculeindex)
    atomicnum=moleculeatom.GetAtomicNum()
    val=moleculeatom.GetExplicitValence()
    atomisinring=moleculeatom.IsInRing()
    probeatom=mol.GetAtom(probeindex)
    moleculeatomcoords=np.array([moleculeatom.GetX(),moleculeatom.GetY(),moleculeatom.GetZ()])
    probeatomcoords=np.array([probeatom.GetX(),probeatom.GetY(),probeatom.GetZ()])
    molprobedist=np.linalg.norm(probeatomcoords-moleculeatomcoords)
    singlepairangle=22.5*0.0174533 # +/- angle from middle
    doublepairangle=60*0.0174533 # +/- angle from middle
    triplepairangle=22.5*0.0174533 # +/- angle from midle pair, just treat as 3 seperate cones
    if poltype.addlonepairvdwsites==True:
        check=False
        if atomicnum==7 and val==3 and atomisinring==False:
            angle=singlepairangle
            atomiter=openbabel.OBAtomAtomIter(moleculeatom)
            natoms=[]
            for natom in atomiter:
                natoms.append(natom)
            npos=[np.array([natom.GetX(),natom.GetY(),natom.GetZ()]) for natom in natoms]
            displacements=[pos-moleculeatomcoords for pos in npos]
            dists=[np.linalg.norm(natomcoords-moleculeatomcoords) for natomcoords in npos]
            lonepairvec=np.zeros(3)
            for i in range(len(displacements)):
                disp=displacements[i]
                dist=dists[i]
                unitvec=disp/dist
                lonepairvec+=unitvec 
                norm=np.linalg.norm(lonepairvec)
                lonepairvec=lonepairvec/norm
            lonepairvec=-1*lonepairvec 
            check=True
            lonepairvecs=[lonepairvec]

        elif atomicnum==7 and val==3 and atomisinring==True:
            angle=singlepairangle
            check=True
            lonepairvecs=LonePairAboveBelowPlane(moleculeatom,mol,moleculeatomcoords,moleculeindex)

        elif atomicnum==7 and val==2 and atomisinring==True:
            angle=singlepairangle
            check=True
            lonepairvecs=LonePairAboveBelowPlane(moleculeatom,mol,moleculeatomcoords,moleculeindex)


        elif atomicnum==8 or atomicnum==16 and atomisinring==True:
            angle=singlepairangle
            check=True
            lonepairvecs=LonePairAboveBelowPlane(moleculeatom,mol,moleculeatomcoords,moleculeindex)


        elif atomicnum==8 or atomicnum==16 and atomisinring==False:
            angle=doublepairangle
            check=True
            lonepairvecs=LonePairBisector(moleculeatom,mol,moleculeatomcoords,moleculeindex)

        elif atomicnum==9 or atomicnum==17 or atomicnum==35 or atomicnum==53:
            angle=triplepairangle
            check=True
            atomiter=openbabel.OBAtomAtomIter(moleculeatom)
            natoms=[]
            for natom in atomiter:
                natoms.append(natom)
            npos=[np.array([natom.GetX(),natom.GetY(),natom.GetZ()]) for natom in natoms]
            displacements=[pos-moleculeatomcoords for pos in npos]
            dists=[np.linalg.norm(natomcoords-moleculeatomcoords) for natomcoords in npos]
            unitvector=displacements[0]/dists[0]
            lonepairvec=-1*unitvector
            ux=lonepairvec[0]
            uy=lonepairvec[1]
            uz=lonepairvec[2]
            x=ux
            y=uy
            z=-(ux*x+uy*y)/uz
            anothervec=np.array([x,y,z])
            norm=np.linalg.norm(anothervec)
            anothervec=anothervec/norm 
            lonepairvecs=[lonepairvec,anothervec]

        if check==True:
            foundone=False
            for i in range(len(lonepairvecs)):
                lonepairvec=lonepairvecs[i]
                x1=moleculeatomcoords
                x2=x1+lonepairvec
                x0=probeatomcoords
                d=np.linalg.norm(np.cross(x2-x1,x1-x0))/np.linalg.norm(x2-x1)
                t=-1*np.dot(x1-x0,x2-x1)/np.linalg.norm(x2-x1)**2
                hvec=x1+t*lonepairvec
                h=np.linalg.norm(hvec-x1)
                R=h*np.tan(angle)
                if R<=d:
                    added=True
                    newdimerfiles.append([newfilename])
                    newprobeindices.append([probeindex])
                    newmoleculeindices.append([moleculeindex])
                    newnumberprobeatoms.append([probeatoms]) 
                if added==True:
                    break


    return newdimerfiles,newprobeindices,newmoleculeindices,newnumberprobeatoms,added



def FindMinimumPoints(poltype,dimerfiles,probeindices,moleculeindices,numberprobeatoms):
    newdimerfiles=[]
    newprobeindices=[]
    newmoleculeindices=[]
    newnumberprobeatoms=[]
    dimertypetofilenamearray={}
    dimertypetoenergyarray={}
    dimertypetoprobeindexarray={}
    dimertypetonumberofprobeatoms={}
    obConversion = openbabel.OBConversion()
    moleculeindextoaddedlonepair={}
    for i in range(len(probeindices)):
        filenameslistoflist=dimerfiles[i]
        probeindexlistoflist=probeindices[i]
        moleculeindexlistoflist=moleculeindices[i]
        probeatoms=numberprobeatoms[i]
        filenameslist = [item for sublist in filenameslistoflist for item in sublist]
        probeindexlist = [item for sublist in probeindexlistoflist for item in sublist]
        moleculeindexlist = [item for sublist in moleculeindexlistoflist for item in sublist]
        for j in range(len(filenameslist)):
            filename=filenameslist[j]
            probeindex=probeindexlist[j]
            moleculeindex=moleculeindexlist[j]
            added=False 
            newfilename=filename.replace('.xyz_2','cart.xyz')
            ConvertTinktoXYZ(poltype,filename,newfilename)
            mol = openbabel.OBMol()
            inFormat = obConversion.FormatFromExt(newfilename)
            obConversion.SetInFormat(inFormat)
            obConversion.ReadFile(mol, newfilename)
            checktoofar=CheckIfProbeIsTooFarOrTooClose(poltype,mol,probeindex,moleculeindex,probeatoms)
            if checktoofar==True:
                continue
            filenameout=filename.replace('.xyz_2','.alz')
            term,error=poltype.CheckNormalTermination(filenameout)
            if term==True and error==False:
                pass
            else:
                cmdstr=poltype.analyzeexe+' '+filename+' '+'-k '+poltype.key4fname +' e '+'> '+filenameout
                poltype.call_subsystem([cmdstr],True)
            energy=ReadIntermolecularEnergyMM(poltype,filenameout)
            if energy==None:
                continue
            hetero=True
            if 'water'not in newfilename:
                hetero=False
            if hetero==True:
                dimertype=tuple([moleculeindex,probeindex,hetero])
            else:
                dimertype=tuple([moleculeindex]) # keep only one of water orientations
            if moleculeindex not in moleculeindextoaddedlonepair.keys():
                moleculeindextoaddedlonepair[moleculeindex]=False
            if moleculeindextoaddedlonepair[moleculeindex]==False:
                newdimerfiles,newprobeindices,newmoleculeindices,newnumberprobeatoms,added=AddLonePairPoints(poltype,newdimerfiles,newprobeindices,newmoleculeindices,newnumberprobeatoms,moleculeindex,probeindex,mol,probeatoms,newfilename)
            if added==False:
                if dimertype not in dimertypetofilenamearray.keys():
                    dimertypetofilenamearray[dimertype]=[]
                dimertypetofilenamearray[dimertype].append(newfilename)
                if dimertype not in dimertypetoenergyarray.keys():
                    dimertypetoenergyarray[dimertype]=[]
                dimertypetoenergyarray[dimertype].append(energy)
                if dimertype not in dimertypetoprobeindexarray.keys():
                    dimertypetoprobeindexarray[dimertype]=[]
                dimertypetoprobeindexarray[dimertype].append(probeindex)
                if dimertype not in dimertypetonumberofprobeatoms.keys():
                    dimertypetonumberofprobeatoms[dimertype]=[]
                dimertypetonumberofprobeatoms[dimertype].append(probeatoms)
            else:
                moleculeindextoaddedlonepair[moleculeindex]=True
    newdimerfiles,newprobeindices,newmoleculeindices,newnumberprobeatoms=SortStructuresViaEnergy(poltype,newdimerfiles,newprobeindices,newmoleculeindices,newnumberprobeatoms,dimertypetoenergyarray,dimertypetofilenamearray,dimertypetoprobeindexarray,dimertypetonumberofprobeatoms)
    return newdimerfiles,newprobeindices,newmoleculeindices,newnumberprobeatoms

def SortStructuresViaEnergy(poltype,newdimerfiles,newprobeindices,newmoleculeindices,newnumberprobeatoms,dimertypetoenergyarray,dimertypetofilenamearray,dimertypetoprobeindexarray,dimertypetonumberofprobeatoms):
    for dimertype,energyarray in dimertypetoenergyarray.items():
        moleculeindex=dimertype[0]
        s = np.array(energyarray)
        sort_index = np.argsort(s)
        filenamearray=np.array(dimertypetofilenamearray[dimertype])
        filenamearray=filenamearray[sort_index]
        probeindexarray=dimertypetoprobeindexarray[dimertype]
        probeindexarray=np.array(probeindexarray)
        probeindexarray=probeindexarray[sort_index]
        probeatomarray=dimertypetonumberofprobeatoms[dimertype]
        probeatomarray=np.array(probeatomarray)
        probeatomarray=probeatomarray[sort_index]
        count=0
        tempnewmoleculeindices=[]
        tempnewprobeindices=[]
        tempnewdimerfiles=[]
        tempprobeatoms=[]
        for i in range(len(filenamearray)):
            filename=filenamearray[i]
            
            probeindex=probeindexarray[i]
            probeatoms=probeatomarray[i]
            count+=1
            if count>poltype.vdwmaxqmstartingpointspertype:
                break

            tempnewmoleculeindices.append(moleculeindex)
            tempnewprobeindices.append(probeindex)
            tempnewdimerfiles.append(filename)
            tempprobeatoms.append(probeatoms)
        newdimerfiles.append(tempnewdimerfiles)
        newprobeindices.append(tempnewprobeindices)
        newmoleculeindices.append(tempnewmoleculeindices)
        newnumberprobeatoms.append(tempprobeatoms)

    return newdimerfiles,newprobeindices,newmoleculeindices,newnumberprobeatoms


def WriteFittingResults(poltype,keyname,classkeytofitresults):
    temp=open(keyname,'r')
    results=temp.readlines()
    temp.close()
    tempname=keyname.replace('.key','_TEMP.key')
    temp=open(tempname,'w')
    for line in results:
        linesplit=line.split()
        if 'vdw ' in line and '#' not in line:
            fwd='%d' % (int(linesplit[1]))       
            for classkey,fitresults in classkeytofitresults.items():
                if classkey==fwd:
                    newresults='# '+fitresults
                    temp.write(newresults)     

        temp.write(line)
    temp.close()
    os.remove(keyname)
    os.rename(tempname,keyname)



def VanDerWaalsOptimization(poltype,missingvdwatomindices):
    poltype.parentdir=os.getcwd()+r'/'
    vdwfoldername='vdw'
    if not os.path.isdir(vdwfoldername):
        os.mkdir(vdwfoldername)
    shutil.copy(poltype.key4fname,vdwfoldername+r'/'+poltype.key4fname)
    CheckAllVdwTypesExist(poltype,poltype.key4fname)
    shutil.copy(poltype.xyzoutfile,vdwfoldername+r'/'+poltype.xyzoutfile)
    shutil.copy(poltype.xyzfname,vdwfoldername+r'/'+poltype.xyzfname)
    os.chdir(vdwfoldername) 
    poltype.optmaxcycle=400
    poltype.optmethod='wB97X-D'
    if poltype.accuratevdwsp==True:
        poltype.espmethod='MP2'
        poltype.espbasisset="aug-cc-pV[TQ]Z"

    else:
        poltype.espmethod='wB97X-D'
        poltype.espbasisset="aug-cc-pVDZ"
    if poltype.debugmode==True:
        poltype.espmethod='hf'
        poltype.espbasisset="minix"

    tempuse_gaus=poltype.use_gaus
    tempuse_gausoptonly=poltype.use_gausoptonly
    if poltype.use_gau_vdw==True and poltype.foundgauss==True:
         poltype.use_gaus=True
    poltype.SanitizeAllQMMethods()
    paramhead=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ParameterFiles/amoebabio18.prm"
    ReplaceParameterFileHeader(poltype,paramhead,poltype.key4fname)
    array=[.8,.85,.9,1,1.1,1.2,1.4]


    dimerfiles,probeindices,moleculeindices,numberprobeatoms=GenerateInitialProbeStructures(poltype,missingvdwatomindices)
    dimerfiles,probeindices,moleculeindices,numberprobeatoms=FindMinimumPoints(poltype,dimerfiles,probeindices,moleculeindices,numberprobeatoms)

    obConversion = openbabel.OBConversion()
    checkarray=[]
    fullprefixarrays=[]
    fulldistarrays=[]
    alloutputfilenames=[]
    originalcharge=poltype.mol.GetTotalCharge()
    originalmul=poltype.mol.GetTotalSpinMultiplicity()
    moleculeprobeindicestoignore=[] # if opt fails
    for i in range(len(probeindices)):
        filenameslist=dimerfiles[i]
        probeindexlist=probeindices[i]
        moleculeindexlist=moleculeindices[i]
        probeatomslist=numberprobeatoms[i]
        distancearrays=[]
        prefixarrays=[]
        filenamesarray=[] 
        for probeidx in range(len(probeindexlist)):
            probeatoms=int(probeatomslist[probeidx])
            filename=filenameslist[probeidx]
            mol = openbabel.OBMol()
            if 'water' in filename:
                totchg=originalcharge
                totmul=originalmul
            else: # homodimer
                totchg=2*originalcharge
                totmul=originalmul# assume is always 1
            mol.SetTotalCharge(totchg)
            mol.SetTotalSpinMultiplicity(totmul)
            chk=mol.GetTotalCharge()
            probeindex=int(probeindexlist[probeidx])
    
            moleculeindex=moleculeindexlist[probeidx]
            inFormat = obConversion.FormatFromExt(filename)
            obConversion.SetInFormat(inFormat)
            obConversion.ReadFile(mol, filename)
            prefix=filename.replace('.xyz','')
            check=CheckIfFittingCompleted(poltype,prefix)
            checkarray.append(check)
            if poltype.use_qmopt_vdw==True:
                poltype.comoptfname=prefix+'-opt.com'
                poltype.chkoptfname=prefix+'-opt.chk'
                poltype.fckoptfname=prefix+'-opt.fchk'
                poltype.logoptfname=prefix+'-opt.log'
                poltype.gausoptfname=prefix+'-opt.log'

                optmol,error = opt.GeometryOptimization(poltype,mol,suffix='1',loose=True,checkbonds=False,modred=False,bondanglerestraints=restraints,skipscferror=False,charge=totchg,skiperrors=True)
                if error==True:
                    moleculeprobeindicestoignore.append([moleculeindex,probeindex])
                    continue
            prefixarrays.append(prefix)
            dimeratoms=mol.NumAtoms()
            moleculeatoms=dimeratoms-probeatoms
            if poltype.use_qmopt_vdw==True:
                moleculeatom=optmol.GetAtom(moleculeindex)
                probeatom=optmol.GetAtom(probeindex)
            else:
                moleculeatom=mol.GetAtom(moleculeindex)
                probeatom=mol.GetAtom(probeindex)

            try:
                moleculeatomcoords=np.array([moleculeatom.GetX(),moleculeatom.GetY(),moleculeatom.GetZ()])
                probeatomcoords=np.array([probeatom.GetX(),probeatom.GetY(),probeatom.GetZ()])
            except:
                moleculeprobeindicestoignore.append([moleculeindex,probeindex])
                continue
            moleculeatomcoords=np.array([moleculeatom.GetX(),moleculeatom.GetY(),moleculeatom.GetZ()])
            probeatomcoords=np.array([probeatom.GetX(),probeatom.GetY(),probeatom.GetZ()])
            equildistance=np.linalg.norm(probeatomcoords-moleculeatomcoords)
            distarray=np.multiply(equildistance,np.array(array))
            distancearrays.append(distarray)
            outputprefixname=filename.split('.')[0]
            outputxyz=outputprefixname+'_tinker.xyz'
            inputxyz=outputprefixname+'_cartesian.xyz'
            if poltype.use_qmopt_vdw==True:
                WriteOutCartesianXYZ(poltype,optmol,inputxyz)
            else:
                WriteOutCartesianXYZ(poltype,mol,inputxyz)

            if 'water' in outputxyz:
                waterbool=True
            else:
                waterbool=False
            ConvertProbeDimerXYZToTinkerXYZ(poltype,inputxyz,poltype.xyzoutfile,outputxyz,waterbool,probeatoms)
            filenamearray=MoveDimerAboutMinima(poltype,outputxyz,outputprefixname,moleculeatoms,moleculeindex,probeindex,equildistance,array)
            qmfilenamearray=GenerateSPInputFiles(poltype,filenamearray,poltype.mol,probeatoms)
            outputfilenames=ExecuteSPJobs(poltype,qmfilenamearray,prefix)
            filenamesarray.append(outputfilenames)
        fullprefixarrays.append(prefixarrays)
        fulldistarrays.append(distancearrays)
        alloutputfilenames.append(filenamesarray)
    dothefit=False
    for check in checkarray:
        if check==False:
            dothefit=True
    classkeytofitresults={}
    if dothefit==True:
        probeindices,moleculeindices=RemoveIgnoredIndices(poltype,probeindices,moleculeindices,moleculeprobeindicestoignore)
        newprobeindices,newmoleculeindices,newprefixarrays,newdistarrays,newoutputfilenames=CombineProbesThatNeedToBeFitTogether(poltype,probeindices,moleculeindices,fullprefixarrays,fulldistarrays,alloutputfilenames)
        for k in range(len(newprobeindices)):
            goodfit=False
            count=1
            subprobeindices=newprobeindices[k]
            submoleculeindices=newmoleculeindices[k]
            subprefixarrays=newprefixarrays[k]
            subdistarrays=newdistarrays[k]
            subfilenames=newoutputfilenames[k]
            flat_probeindices = [item for sublist in subprobeindices for item in sublist]
            flat_moleculeindices = [item for sublist in submoleculeindices for item in sublist]
            flat_prefixarrays = [item for sublist in subprefixarrays for item in sublist]
            flat_distarrays = [item for sublist in subdistarrays for item in sublist]
            flat_filenames = [item for sublist in subfilenames for item in sublist]
            flat_filenames = [item for sublist in flat_filenames for item in sublist]
            while goodfit==False:
                if count>2:
                    break
                vdwtypesarray=[]
                initialradii=[]
                initialdepths=[]
                initialreds=[]
                minradii=[]
                maxradii=[]
                mindepths=[]
                maxdepths=[]
                minreds=[]
                maxreds=[]
                fitredboolarray=[]
                ReadCounterPoiseAndWriteQMData(poltype,flat_filenames)
                for probeidx in range(len(flat_probeindices)):
                    prefix=flat_prefixarrays[probeidx]
                    probeindex=flat_probeindices[probeidx]
                    moleculeindex=flat_moleculeindices[probeidx]
                    vdwtype=poltype.idxtosymclass[moleculeindex]
                    atom=poltype.mol.GetAtom(moleculeindex)
                    valence=atom.GetExplicitValence()
                    atomicnum=atom.GetAtomicNum()
                    if valence==1 and atomicnum!=8 and atomicnum!=16 and poltype.fitred==True:
                       fitred=True
                    else:
                       fitred=False
                    if vdwtype not in vdwtypesarray:
                        vdwtypesarray.append(vdwtype)
                        fitredboolarray.append(fitred)
                        initialvdwradius,initialvdwdepth,minvdwradius,maxvdwradius,minvdwdepth,maxvdwdepth,red,minred,maxred=GrabVdwParameters(poltype,vdwtype)
                        initialradii.append(initialvdwradius)
                        initialdepths.append(initialvdwdepth)
                        minradii.append(minvdwradius)
                        maxradii.append(maxvdwradius)
                        mindepths.append(minvdwdepth)
                        maxdepths.append(maxvdwdepth)
                        initialreds.append(red)
                        minreds.append(minred)
                        maxreds.append(maxred)
                    if 'water' not in prefix:
                        adjustedprobeindex=int(probeindex-len(poltype.idxtosymclass.keys()))
                        vdwtype=poltype.idxtosymclass[adjustedprobeindex]
                        atom=poltype.mol.GetAtom(adjustedprobeindex)
                        valence=atom.GetExplicitValence()
                        if valence==1:
                           fitred=True
                        else:
                           fitred=False

                        if vdwtype not in vdwtypesarray and adjustedprobeindex in missingvdwatomindices:
                            vdwtypesarray.append(vdwtype)
                            fitredboolarray.append(fitred)
                            initialvdwradius,initialvdwdepth,minvdwradius,maxvdwradius,minvdwdepth,maxvdwdepth,red,minred,maxred=GrabVdwParameters(poltype,vdwtype)
                            initialradii.append(initialvdwradius)
                            initialdepths.append(initialvdwdepth)
                            minradii.append(minvdwradius)
                            maxradii.append(maxvdwradius)
                            mindepths.append(minvdwdepth)
                            maxdepths.append(maxvdwdepth)
                            initialreds.append(red)
                            minreds.append(minred)
                            maxreds.append(maxred)
                WriteInitialPrmFile(poltype,vdwtypesarray,initialradii,initialdepths,minradii,maxradii,mindepths,maxdepths,initialreds,minreds,maxreds)
                vdwradii,vdwdepths,vdwreds=VDWOptimizer(poltype,count,fitredboolarray)
                for k in range(len(flat_prefixarrays)):
                    prefix=flat_prefixarrays[k]
                    distarray=flat_distarrays[k]
                    PlotEnergyVsDistance(poltype,distarray,prefix,vdwradii,vdwdepths,vdwreds,vdwtypesarray,count)
                    othergoodfit,classkeytofitresults=PlotQMVsMMEnergy(poltype,vdwtypesarray,prefix,count,classkeytofitresults)
                goodfit,classkeytofitresults=PlotQMVsMMEnergy(poltype,vdwtypesarray,flat_prefixarrays,count,classkeytofitresults,allprefix=True)
                count+=1
    WriteFittingResults(poltype,poltype.key4fname,classkeytofitresults)
    shutil.copy(poltype.key4fname,'../'+poltype.key5fname)
    os.chdir(poltype.parentdir)
    poltype.use_gaus=tempuse_gaus
    poltype.use_gausoptonly=tempuse_gausoptonly

