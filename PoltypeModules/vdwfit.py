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
import shutil
import copy

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
    return 1 - (squared_error_regr/squared_error_y_mean)


def MeanError(data,pred):
    Sum=0
    for i in range(len(data)):
        true=data[i]
        predv=pred[i]
        diff=true-predv
        Sum+=diff
    Sum=Sum/len(data)
    return Sum
    

def ReadInitialPrmFile(poltype):
    temp=open('INITIAL.PRM','r')
    frontarray=[]
    for line in temp.readlines():
        linesplit=re.split(r'(\s+)', line)
        front=linesplit[0:3]
        frontarray.append(''.join(front))
    temp.close()
    return frontarray


def RMSE(rawdata,refdata):
    rms=0.0
    len1=len(rawdata)
    len2=len(refdata)
    if len1 != len2:
        print("Data count does not match!")
    else:
        for n in range(len1):
          rms=rms+(float(rawdata[n])-float(refdata[n]))**2.0
        rms=(rms/len1)**0.5
    return rms

def writePRM(poltype,params,vdwtype):
    for i in range(0,len(params),2):
        oFile = open("temp.key", 'w')
        rFile=open(poltype.key5fname,'r')
        for line in rFile.readlines():
            if vdwtype in line and 'vdw' in line:
                oFile.write('vdw '+vdwtype+" %s %s\n"%(params[i], params[i+1]))
            else:
                oFile.write(line)
    oFile.flush()
    os.fsync(oFile.fileno())
    oFile.close()
    rFile.close()
    os.remove(poltype.key5fname)
    os.rename("temp.key",poltype.key5fname)
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
     

def myFUNC(params,poltype,vdwtype):
    params=list(params)
    writePRM(poltype,params,vdwtype)
    target = []
    target = readOneColumn("QM_DATA",1)
    target=[float(i) for i in target]
    target=[i-min(target) for i in target]
    temp=open('QM_DATA','r')
    cmdarray=[] 
    filenamearray=[]
    for line in temp.readlines():
        xyzname=line.split()[0]
        filename=xyzname.replace('.xyz','.alz')
        cmdstr=poltype.analyzeexe+' '+xyzname+' '+'-k '+poltype.key5fname +' e '+'> '+filename
        cmdarray.append(cmdstr)
        filenamearray.append(filename)

    temp.close()
    # need to execute commands now
    for cmd in cmdarray:
        poltype.call_subsystem(cmd,True)

    ReadAnalyzeEnergiesWriteOut(poltype,filenamearray)
    vdw=readOneColumn("SP.dat",-1)
    current=list(np.array(vdw))
    current=[float(i) for i in current]
    current=[i-min(current) for i in current]
    new_rmse=RMSE(current,target)

    return new_rmse


def PlotQMVsMMEnergy(poltype,vdwtypesarray,prefix):
    target = readOneColumn("QM_DATA",1,prefix)
    target=[float(i) for i in target]
    target=[i-min(target) for i in target]
    vdw=readOneColumn("SP.dat",-1,prefix)
    current=list(np.array(vdw))
    current=[float(i) for i in current]
    current=[i-min(current) for i in current]
    vdwtypes=[str(i) for i in vdwtypesarray]
    vdwtypestring=','.join(vdwtypes)

    MSE=MeanError(current,target)
    MAE=metrics.mean_absolute_error(current,target)
    m, b = best_fit_slope_and_intercept(current,target)
    regression_line = [(m*x)+b for x in current]
    new_rmse=RMSE(current,target)

    r_squared = coefficient_of_determination(current,target)
    fig = plt.figure()
    plt.plot(current,target,label='R^2=%s MAE=%s RMSE=%s MSE=%s'%(round(r_squared,2),round(MAE,2),round(new_rmse,2),round(MSE,2)))
    plt.plot(current,regression_line,label='Linear Regression line')
    plt.ylabel('QM BSSE Corrected (kcal/mol)')
    plt.xlabel('AMOEBA (kcal/mol)')
    plt.legend(loc='best')
    plt.title('QM vs AMOEBA , '+vdwtypestring)
    fig.savefig('QMvsAMOEBA-'+prefix+'_'+vdwtypestring+'.png')


def VDWOptimizer(poltype):
    x0 = []
    temp=open("INITIAL.PRM",'r')
    lines = temp.readlines()
    for line in lines:
        x0.append(float(line.split()[2]))
        x0.append(float(line.split()[3]))
    x0 = np.array(x0)
    rmax=readOneColumn("INITIAL.PRM", 5)
    rmin=readOneColumn("INITIAL.PRM", 4)
    depthmax=readOneColumn("INITIAL.PRM", 7)
    depthmin=readOneColumn("INITIAL.PRM", 6)
    vdwtypes=readOneColumn("INITIAL.PRM", 1)
    vdwtype=vdwtypes[0]
    l1=list(zip(rmin, rmax))
    l2=list(zip(depthmin, depthmax))
    MyBounds=[]
    for i in range(len(l1)):
        
        l1bound=l1[i]
        l2bound=l2[i]
        MyBounds.append(l1bound)
        MyBounds.append(l2bound)
    tuple(MyBounds)
    ''' local optimization method can be BFGS, CG, Newton-CG, L-BFGS-B,etc.., see here\
    https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.optimize.minimize.html'''
    errorfunc= lambda p: (myFUNC(p,poltype,vdwtype))

    ret = minimize(myFUNC, x0, method='L-BFGS-B', jac=None, bounds=MyBounds, args=(poltype,vdwtype),options={'gtol': 1e-1, 'eps':1e-4, 'disp': True})
    vdwradius=round(ret.x[0],3)
    vdwdepth=round(ret.x[1],3)
    ofile = open("RESULTS.OPT", "a")
    for i in range(0,len(ret.x),2):
        ofile.write("Param:%15.6f%15.6f\n"%(ret.x[i], ret.x[i+1]))
    ofile.write("%s\n"%(ret.message))
    ofile.write("\n")
    ofile.close()
    return vdwradius,vdwdepth 

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


def GrabMonomerEnergy(poltype,line):
    linesplit=line.split()
    energy=float(linesplit[4]) 
    monomerenergy=(energy)*poltype.Hartree2kcal_mol
    return monomerenergy



def ReadCounterPoiseAndWriteQMData(poltype,logfilelist):
    temp=open('QM_DATA','w')
    for f in logfilelist:
        tmpfh=open(f,'r')
        if poltype.use_gaus==True:
            frag1calc=False
            frag2calc=False
            for line in tmpfh:
                if 'Counterpoise: corrected energy =' in line:
                    dimerenergy=float(line.split()[4])*poltype.Hartree2kcal_mol
                elif 'Counterpoise: doing DCBS calculation for fragment   1' in line:
                    frag1calc=True
                    frag2calc=False
                elif 'Counterpoise: doing DCBS calculation for fragment   2' in line:
                    frag1calc=False
                    frag2calc=True
                elif 'SCF Done' in line:
                    if frag1calc:
                        frag1energy=GrabMonomerEnergy(poltype,line)
                    elif frag2calc:
                        frag2energy=GrabMonomerEnergy(poltype,line)
            interenergy=dimerenergy-(frag1energy+frag2energy)
        else:
            for line in tmpfh:
                if 'CP Energy =' in line and 'print' not in line:
                    linesplit=line.split()
                    interenergy=float(linesplit[3])*poltype.Hartree2kcal_mol

        tmpfh.close()
        temp.write(f.replace('_sp.log','.xyz')+' '+str(interenergy)+'\n')
    temp.close()




def PlotEnergyVsDistance(poltype,distarray,prefix,rad,depth,vdwtypesarray):
    vdwtypes=[str(i) for i in vdwtypesarray]
    vdwtypestring=','.join(vdwtypes)
    qmenergyarray = readOneColumn("QM_DATA",1,prefix)
    qmenergyarray=[float(i) for i in qmenergyarray]
    qmenergyarray=[i-min(qmenergyarray) for i in qmenergyarray]
    vdw=readOneColumn("SP.dat",-1,prefix)
    energyarray=list(np.array(vdw))
    energyarray=[float(i) for i in energyarray]
    energyarray=[i-min(energyarray) for i in energyarray]
    def RMSD(c):
        return np.sqrt(np.mean(np.square(np.add(np.subtract(np.array(energyarray),np.array(qmenergyarray)),c))))


    r_squared = round(coefficient_of_determination(energyarray,qmenergyarray),2)
    result=fmin(RMSD,.5)
    minRMSD=round(RMSD(result[0]),2)
    plotname='EnergyVsDistance-'+prefix+'_'+vdwtypestring+'.png'
    fig = plt.figure()
    title=prefix+' VdwTypes = '+vdwtypestring
    plt.title(title)
    plt.plot(distarray,energyarray,'b-',label='MM ,'+'Radius=%s, Depth=%s'%(round(rad,2),round(depth,2)))
    plt.plot(distarray,qmenergyarray,'r-',label='QM')
    plt.plot()
    plt.ylabel('Energy (kcal/mol)')
    plt.xlabel('Distance Angstrom '+'RMSD=%s, R^2=%s'%(minRMSD,r_squared))
    plt.legend(loc='best')
    fig.savefig(plotname)


def ReadIntermolecularEnergyMM(poltype,filename):
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


def WriteInitialPrmFile(poltype,vdwtypesarray,initialradii,initialdepths,minradii,maxradii,mindepths,maxdepths):
    temp=open('INITIAL.PRM','w')
    for i in range(len(vdwtypesarray)):
        vdwtype=str(vdwtypesarray[i])
        initialradius=str(initialradii[i])
        initialdepth=str(initialdepths[i])
        minradius=str(minradii[i])
        maxradius=str(maxradii[i])
        mindepth=str(mindepths[i])
        maxdepth=str(maxdepths[i])
        line='vdw '+vdwtype+' '+initialradius+' '+initialdepth+' '+minradius+' '+maxradius+' '+mindepth+' '+maxdepth+'\n'
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

    for n in range(len(atoms)):
        if n>=len(atoms)-probeatoms:
            tmpfh.write("%3s%s             %14.7f%14.7f%14.7f\n"%(atoms[n],'(Fragment=2)',float(coord[n][0]),float(coord[n][1]),float(coord[n][2]))) 
        else:
            tmpfh.write("%3s%s             %14.7f%14.7f%14.7f\n"%(atoms[n],'(Fragment=1)',float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))    
    
    if ('I ' in poltype.mol.GetSpacedFormula()):
        formulalist=poltype.mol.GetSpacedFormula().lstrip().rstrip().split()
        temp=open(poltype.basissetpath+basissetfile,'r')
        results=temp.readlines()
        temp.close()
        for line in results:
            if '!' not in line:
                tmpfh.write(line)

        temp=open(poltype.basissetpath+iodinebasissetfile,'r')
        results=temp.readlines()
        temp.close()
        for line in results:
            if '!' not in line:
                tmpfh.write(line)


        temp=open(poltype.basissetpath+iodinebasissetfile,'r')
        results=temp.readlines()
        temp.close()
        for line in results:
            if '!' not in line:
                tmpfh.write(line)


    tmpfh.write("\n")
    tmpfh.close()



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
    spacedformulastr=mol.GetSpacedFormula()
    if ('I ' in spacedformulastr):
        temp.write('basis {'+'\n')
        temp.write('['+' '+poltype.espbasissetfile+' '+poltype.iodineespbasissetfile +' '+ ']'+'\n')
        temp=ReadInBasisSet(poltype,temp,poltype.espbasissetfile,poltype.iodineespbasissetfile)
        temp.write('}'+'\n')
        temp.write("e_dim= energy('%s',bsse_type='cp')" % (poltype.espmethod.lower())+'\n')
    else:
        temp.write("e_dim= energy('%s/%s',bsse_type='cp')" % (poltype.espmethod.lower(),poltype.espbasisset)+'\n')
    temp.write('\n')
    temp.write('molecule { '+'\n')
    temp.write('%d %d\n' % (chg, 1))
    for n in range(len(atoms)):
        if n>=len(atoms)-probeatoms:
            temp.write("%3s             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2]))) 
        else:
            temp.write("%3s             %14.7f%14.7f%14.7f\n"%('@'+atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))    
    temp.write('}'+'\n')

    if ('I ' in spacedformulastr):
        temp.write('basis {'+'\n')
        temp.write('['+' '+poltype.espbasissetfile+' '+poltype.iodineespbasissetfile +' '+ ']'+'\n')
        temp=ReadInBasisSet(poltype,temp,poltype.espbasissetfile,poltype.iodineespbasissetfile)
        temp.write('}'+'\n')
        temp.write("e_mon_a= energy('%s')" % (poltype.espmethod.lower())+'\n')
    else:
        temp.write("e_mon_a= energy('%s/%s')" % (poltype.espmethod.lower(),poltype.espbasisset)+'\n')

    temp.write('\n')
    temp.write('molecule { '+'\n')
    temp.write('%d %d\n' % (chg, 1))
    for n in range(len(atoms)):
        if n>=len(atoms)-probeatoms:
            temp.write("%s             %14.7f%14.7f%14.7f\n"%('@'+atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2]))) 
        else:
            temp.write("%3s             %14.7f%14.7f%14.7f\n"%(atoms[n],float(coord[n][0]),float(coord[n][1]),float(coord[n][2])))    
    temp.write('}'+'\n')

    if ('I ' in spacedformulastr):
        temp.write('basis {'+'\n')
        temp.write('['+' '+poltype.espbasissetfile+' '+poltype.iodineespbasissetfile +' '+ ']'+'\n')
        temp=ReadInBasisSet(poltype,temp,poltype.espbasissetfile,poltype.iodineespbasissetfile)
        temp.write('}'+'\n')
        temp.write("e_mon_b= energy('%s')" % (poltype.espmethod.lower())+'\n')
    else:
        temp.write("e_mon_b= energy('%s/%s')" % (poltype.espmethod.lower(),poltype.espbasisset)+'\n')

    temp.write('\n')
    
    temp.write("e_cp = e_dim - e_mon_a - e_mon_b"+'\n')
    temp.write("psi4.print_out('CP Energy = %10.6f' % (e_cp))"+'\n')
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
    etab = openbabel.OBElementTable()
    for atom in atomiter:
        atomcounter+=1
        atomsymb=etab.GetSymbol(atom.GetAtomicNum())
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
    outputfilenames=[]
    for i in range(len(qmfilenamearray)):
        filename=qmfilenamearray[i]
        if poltype.use_gaus==True:
            cmdstr = 'cd '+shlex.quote(os.getcwd())+' && '+'GAUSS_SCRDIR='+poltype.scrtmpdirgau+' '+poltype.gausexe+' '+filename
        
            outputname=filename.replace('.com','.log')
        else:
            outputname=filename.replace('.psi4','.log')
            cmdstr='cd '+shlex.quote(os.getcwd())+' && '+'psi4 '+filename+' '+outputname


        
        finished,error=poltype.CheckNormalTermination(outputname)
        if finished==True and error==False:
            pass
        else:
            if os.path.isfile(outputname):
                statinfo=os.stat(outputname)
                size=statinfo.st_size
                if size!=0:
                    os.remove(outputname)
                    listofjobs.append(cmdstr)
                    jobtooutputlog[cmdstr]=os.getcwd()+r'/'+outputname
            else:
                listofjobs.append(cmdstr)
                jobtooutputlog[cmdstr]=os.getcwd()+r'/'+outputname
        outputfilenames.append(outputname)
    lognames=[]
    for job in listofjobs:
        log=jobtooutputlog[job]
        lognames.append(os.path.abspath(poltype.logfname))
    fulljobtolog=dict(zip(listofjobs, lognames)) 
    fulljobtooutputlog.update(jobtooutputlog)

    scratchdir=poltype.scrtmpdirgau
    jobtologlistfilenameprefix=os.getcwd()+r'/'+'QMSPJobToLog'+'_'+prefix
    if poltype.externalapi!=None:
        if len(listofjobs)!=0:
            call.CallExternalAPI(poltype,fulljobtolog,jobtologlistfilenameprefix,scratchdir)
        finshedjobs,errorjobs=poltype.WaitForTermination(fulljobtooutputlog)
    else:
        finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(fulljobtooutputlog,True)

    return outputfilenames


def GrabVdwParameters(poltype,vdwtype):
    temp=open(poltype.key5fname,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'vdw' in line:
            if str(vdwtype) in line:
                linesplit=line.split()
                radius=float(linesplit[2])
                minvdwradius=radius-.1*radius
                maxvdwradius=radius+.1*radius
                depth=float(linesplit[3])
                minvdwdepth=depth-.1*depth
                maxvdwdepth=depth+.1*depth

                return radius,depth,minvdwradius,maxvdwradius,minvdwdepth,maxvdwdepth 


def GenerateCoordinateGuesses(indextoreferencecoordinate):
    coordinatesguess=[]
    for i in range(len(indextoreferencecoordinate.keys())):
        coordinate=indextoreferencecoordinate[i]
        x,y,z=coordinate[:]
        coordinatesguess.append(x)
        coordinatesguess.append(y)
        coordinatesguess.append(z)
    return coordinatesguess

def UpdateCoordinates(coords,indextoreferencecoordinate):
    for i in range(len(indextoreferencecoordinate.keys())):
        startindex=3*i
        coordinate=np.array([coords[startindex],coords[startindex+1],coords[startindex+2]])
        indextoreferencecoordinate[i]=coordinate
    return indextoreferencecoordinate

def GenerateReferenceDistances(indextoreferencecoordinate,indextomolecule,indextotargetatom,indextoreferenceelement,vdwradius):
    indexpairtoreferencedistance={}
    indexpairtobounds={}
    indices=list(indextoreferencecoordinate.keys())
    allpairs=list(itertools.combinations(indices, 2)) 
    pairwiseratio=1
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
                targetdistance=pairwiseratio*(vdwradius1+vdwradius2)
                bound=[targetdistance,targetdistance]

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
                targetdistance=pairwiseratio*(vdwradius1+vdwradius2)
                bound=[targetdistance,targetdistance]
            else:
                bound=[targetdistance,100] # high upper bound, large range to not add cost in cost function                

        indexpairtoreferencedistance[tuple(pair)]=targetdistance 
        indexpairtobounds[tuple(pair)]=bound

 
    return indexpairtoreferencedistance,indexpairtobounds


def GenerateReferenceAngles(poltype,p2,atoms2,p1,atoms1,mol,probemol,indextoreferencecoordinate):
    indicestoreferenceangleprobe={}
    indicestoreferenceanglemoleculeneighb={}
    indicestoreferenceanglemoleculeneighbneighb={}
    acceptorindex=p1
    shiftedp2=len(atoms1)+p2
    donorindex=shiftedp2
    donorneighbs=[]
    probeidxtosymclass,symmetryclass=symm.gen_canonicallabels(poltype,probemol) 
    acceptoratom=mol.GetAtom(p1+1)
    donoratom=probemol.GetAtom(p2+1)
    probeneighbs=[]
    atomatomiter=openbabel.OBAtomAtomIter(donoratom)
    for atom in atomatomiter:
        probeneighbs.append(atom.GetIdx())
    acceptorneighbs=[]
    atomatomiter=openbabel.OBAtomAtomIter(acceptoratom)
    for atom in atomatomiter:
        acceptorneighbs.append(atom.GetIdx())

    probeneighbsymclasses=[probeidxtosymclass[i] for i in probeneighbs] 
    probeneighbsymclassesset=list(set(probeneighbsymclasses))
    probeneighbatoms=[probemol.GetAtom(i) for i in probeneighbs]
    acceptorneighbatoms=[mol.GetAtom(i) for i in acceptorneighbs]

    probeneighbs=[i-1+len(atoms1) for i in probeneighbs] # shift to 0 index, shift passed first molecule
    acceptorneighbs=[i-1 for i in acceptorneighbs]
    donorcoordinate=indextoreferencecoordinate[donorindex]
    acceptorcoordinate=indextoreferencecoordinate[acceptorindex]
    for donorneighbindex in probeneighbs:
        if len(probeneighbs)==2 and len(probeneighbsymclassesset)==1:
            angleatoms=[probeneighbatoms[0],donoratom,probeneighbatoms[1]]
            bisectangle=probemol.GetAngle(angleatoms[0],angleatoms[1],angleatoms[2])
            angle=180-.5*bisectangle 

        elif len(probeneighbs)==1:
            angle=180
        else:

            donorneighbcoordinate=indextoreferencecoordinate[donorneighbindex]
            donortoacceptor=acceptorcoordinate-donorcoordinate
            donortodonorneighb=donorneighbcoordinate-donorcoordinate
            donortoacceptornormed=donortoacceptor/np.linalg.norm(donortoacceptor)
            donortodonorneighbnormed=donortodonorneighb/np.linalg.norm(donortodonorneighb) 
            angle=np.arccos(np.dot(donortoacceptornormed,donortodonorneighbnormed))


        angleindices=tuple([donorneighbindex,donorindex,acceptorindex])
        indicestoreferenceangleprobe[angleindices]=angle
    if len(acceptorneighbs)!=1:
        angle=90
        for acceptorneighbindex in acceptorneighbs:
            angleindices=tuple([donorindex,acceptorneighbindex,acceptorindex])
            indicestoreferenceanglemoleculeneighb[angleindices]=angle
    elif len(acceptorneighbs)==1:
        neighb=acceptorneighbatoms[0]
        atomatomiter=openbabel.OBAtomAtomIter(neighb)
        neighbofneighbsatoms=[]
        for atom in atomatomiter:
            neighbofneighbsatoms.append(atom)
        for atom in neighbofneighbsatoms:
            atomidx=atom.GetIdx()
            angleatoms=[atom,neighb,acceptoratom] 
            angle=mol.GetAngle(angleatoms[0],angleatoms[1],angleatoms[2])
            angleindices=tuple([atomidx-1,neighb.GetIdx()-1,acceptorindex,donorindex]) # acceptorindex and donorindex colinear here, but need both to get correct distance information for law of cosines
            indicestoreferenceanglemoleculeneighbneighb[angleindices]=angle
    
    return indicestoreferenceangleprobe,indicestoreferenceanglemoleculeneighb,indicestoreferenceanglemoleculeneighbneighb


def ConvertAngleRestraintToDistanceRestraint(indexpairtoreferencedistance,indicestoreferenceangleprobe,indicestoreferenceanglemoleculeneighb,indicestoreferenceanglemoleculeneighbneighb,indexpairtobounds,indextoreferencecoordinate):
    for indices,targetangle in indicestoreferenceangleprobe.items(): 
        donorneighbindex=indices[0]
        acceptorindex=indices[2]
        donorindex=indices[1]
        inputpair=tuple([donorneighbindex,donorindex])
        inputdistance=GrabPairwiseDistance(inputpair,indexpairtoreferencedistance)
        targetpair=tuple([donorindex,acceptorindex])
        targetdistance=GrabPairwiseDistance(targetpair,indexpairtoreferencedistance)
        angledist=LawOfCosines(inputdistance,targetdistance,targetangle)    
        anglepair=tuple([donorneighbindex,acceptorindex])
        indexpairtoreferencedistance[anglepair]=angledist
        indexpairtobounds[anglepair]=[angledist,angledist]

    for indices,targetangle in indicestoreferenceanglemoleculeneighb.items(): 
        acceptorneighbindex=indices[1]
        acceptorindex=indices[2]
        donorindex=indices[0]
        inputpair=tuple([acceptorneighbindex,acceptorindex])
        inputdistance=GrabPairwiseDistance(inputpair,indexpairtoreferencedistance)
        targetpair=tuple([donorindex,acceptorindex])
        targetdistance=GrabPairwiseDistance(targetpair,indexpairtoreferencedistance)
        acceptorcoordinate=indextoreferencecoordinate[acceptorindex]
        acceptorneighbcoordinate=indextoreferencecoordinate[acceptorneighbindex]
        donorcoordinate=indextoreferencecoordinate[donorindex]
        donoracceptorvector=acceptorcoordinate-donorcoordinate
        donoracceptorneighbvector=acceptorneighbcoordinate-donorcoordinate
        donoracceptorvectornormed=donoracceptorvector/np.linalg.norm(donoracceptorvector)
        donoracceptorneighbvectornormed=donoracceptorneighbvector/np.linalg.norm(donoracceptorneighbvector)
        currentangle=np.arccos(np.dot(donoracceptorvectornormed,donoracceptorneighbvectornormed))
        angle=180-currentangle-targetangle
        angledist=LawOfCosines(inputdistance,targetdistance,angle)    
        anglepair=tuple([acceptorneighbindex,donorindex])
        indexpairtoreferencedistance[anglepair]=angledist
        indexpairtobounds[anglepair]=[angledist,angledist]

    for indices,targetangle in indicestoreferenceanglemoleculeneighbneighb.items(): 
        acceptorneighbneighbindex=indices[0]
        acceptorneighbindex=indices[1]
        acceptorindex=indices[2]
        donorindex=indices[3] 
        inputpair=tuple([acceptorneighbneighbindex,acceptorneighbindex])
        inputdistance=GrabPairwiseDistance(inputpair,indexpairtoreferencedistance)
        targetpair=tuple([donorindex,acceptorindex])
        targetdistance=GrabPairwiseDistance(targetpair,indexpairtoreferencedistance)
        anotherinputpair=tuple([acceptorneighbindex,acceptorindex])
        anotherinputdistance=GrabPairwiseDistance(anotherinputpair,indexpairtoreferencedistance)
        firstdistance=inputdistance
        seconddistance=targetdistance+anotherinputdistance
        angledist=LawOfCosines(firstdistance,seconddistance,targetangle)    
        anglepair=tuple([acceptorneighbneighbindex,donorindex])
        indexpairtoreferencedistance[anglepair]=angledist
        indexpairtobounds[anglepair]=[angledist,angledist]


           
    return indexpairtoreferencedistance,indexpairtobounds



def GrabPairwiseDistance(pair,indexpairtoreferencedistance):
    if pair in indexpairtoreferencedistance.keys():
        distance=indexpairtoreferencedistance[pair]
    elif pair[::-1] in indexpairtoreferencedistance.keys():
        distance=indexpairtoreferencedistance[pair[::-1]]
    return distance


def LawOfCosines(a,b,angleC):
    return np.sqrt(a**2+b**2-2*a*b*np.cos(np.radians(angleC)))

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


def optimize(poltype,atoms1, atoms2, coords1, coords2, p1, p2, dimer,vdwradius,mol,probemol):
    indextoreferencecoordinate,indextoreferenceelement,indextomolecule,indextotargetatom=GenerateInitialDictionaries(coords1,coords2,atoms1,atoms2,p1,p2) 
    indexpairtoreferencedistance,indexpairtobounds=GenerateReferenceDistances(indextoreferencecoordinate,indextomolecule,indextotargetatom,indextoreferenceelement,vdwradius)
    indicestoreferenceangleprobe,indicestoreferenceanglemoleculeneighb,indicestoreferenceanglemoleculeneighbneighb=GenerateReferenceAngles(poltype,p2,atoms2,p1,atoms1,mol,probemol,indextoreferencecoordinate)
    indexpairtoreferencedistance,indexpairtobounds=ConvertAngleRestraintToDistanceRestraint(indexpairtoreferencedistance,indicestoreferenceangleprobe,indicestoreferenceanglemoleculeneighb,indicestoreferenceanglemoleculeneighbneighb,indexpairtobounds,indextoreferencecoordinate) 
    coordinatesguess=GenerateCoordinateGuesses(indextoreferencecoordinate)
    def PairwiseCostFunction(x):
        func=0
        for indexpair,bounds in indexpairtobounds.items():
            firstindex=indexpair[0]
            secondindex=indexpair[1]
            startfirstindex=3*firstindex
            startsecondindex=3*secondindex     
            firstcoordinate=np.array([x[startfirstindex],x[startfirstindex+1],x[startfirstindex+2]])
            secondcoordinate=np.array([x[startsecondindex],x[startsecondindex+1],x[startsecondindex+2]])       
            distance=np.linalg.norm(firstcoordinate-secondcoordinate)
            referencedistance=indexpairtoreferencedistance[indexpair]
            difference=np.abs(distance-referencedistance)
            lowerbound=bounds[0]
            upperbound=bounds[1]
            if distance<lowerbound or distance>upperbound:
                func+=difference**2


        return func


    sol = minimize(PairwiseCostFunction, coordinatesguess, method='SLSQP',options={'disp':False, 'maxiter': 1000, 'ftol': 1e-6})
    coords=sol.x
    indextoreferencecoordinate=UpdateCoordinates(coords,indextoreferencecoordinate)

    with open(dimer, "w") as f:
      f.write(str(len(indextoreferencecoordinate.keys()))+"\n")
      f.write("\n")
      for index,coordinate in indextoreferencecoordinate.items():
          element=indextoreferenceelement[index]
          f.write("%3s %12.5f%12.5f%12.5f\n"%(element, coordinate[0], coordinate[1], coordinate[2]))


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
        valence=atom.GetValence()
        hyb=atom.GetHyb() 
        if valence==4 and hyb==3:
            indicestodelete.append(i)
    for index in indicestodelete:
        del indexarray[index]

    return indexarray


def GenerateInitialProbeStructure(poltype,missingvdwatomindices):
    molecules = [poltype.xyzfname]
    probes = GenerateProbePathNames(poltype,poltype.vdwprobenames,poltype.xyzfname)
    vdwradius = {"H" : 1.20, "Li": 1.82, "Na": 2.27, "K": 2.75, "Rb": 3.03, "Cs": 3.43, \
                 "Be": 1.53, "Mg": 1.73, "Ca": 2.31, "B": 1.92, "C": 1.70, "N": 1.55, "O":1.52, \
                 "P" : 1.80, "S" : 1.80, "F" : 1.47, "Cl":1.75, "Br":1.85, "Zn":1.39}           
    dimernames=[]
    probeindices=[]
    moleculeindices=[]
    numberprobeatoms=[]
    for mol in molecules:
        atoms1,coords1,order1, types1, connections1=readTXYZ(poltype,mol)
        mol_spots = missingvdwatomindices.copy()
        mol_spots = CheckBuriedAtoms(poltype,mol_spots,poltype.mol)
        mol_spots = [i-1 for i in mol_spots]
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
            prob_spots = CheckBuriedAtoms(poltype,prob_spots,probemol,zeroindex=True)
            for p1 in mol_spots:
                probelist=[]
                probeindiceslist=[]
                moleculeindiceslist=[]
                for p2 in prob_spots:
                    e1=atoms1[p1]
                    e2=atoms2[p2]
                    if e1=='H' and e2=='H':
                        continue
                    dimer = mol[:-4] + "-" + probename[:-4] + "_" + str("%d_%d"%(p1+1,len(atoms1)+p2+1)) + ".xyz"
                    if not os.path.isfile(dimer):
                        optimize(poltype,atoms1, atoms2, coords1, coords2, p1, p2, dimer,vdwradius,poltype.mol,probemol)
                    probelist.append(dimer)
                    probeindiceslist.append(p2+1+len(atoms1))
                    moleculeindiceslist.append(p1+1)
                moleculeindices.append(moleculeindiceslist)
                probeindices.append(probeindiceslist)
                dimernames.append(probelist)
                numberprobeatoms.append(probemol.NumAtoms())
    return dimernames,probeindices,moleculeindices,numberprobeatoms

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
        if 'parameter' in line:
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
 

def VanDerWaalsOptimization(poltype,missingvdwatomindices):
    poltype.parentdir=os.getcwd()+r'/'
    vdwfoldername='vdw'
    if not os.path.isdir(vdwfoldername):
        os.mkdir(vdwfoldername)
    shutil.copy(poltype.key5fname,vdwfoldername+r'/'+poltype.key5fname)
    shutil.copy(poltype.xyzoutfile,vdwfoldername+r'/'+poltype.xyzoutfile)
    shutil.copy(poltype.xyzfname,vdwfoldername+r'/'+poltype.xyzfname)
    os.chdir(vdwfoldername) 
    poltype.optmaxcycle=400
    poltype.optmethod='wB97X-D'
    poltype.espmethod='wB97X-D'
    poltype.espbasisset="aug-cc-pVDZ"
    poltype.use_gaus=False
    poltype.use_gausoptonly=False
    poltype.SanitizeAllQMMethods()
    paramhead=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ParameterFiles/amoebabio18.prm"
    ReplaceParameterFileHeader(poltype,paramhead,poltype.key5fname)
    array=[.8,.9,1,1.1,1.2]
    dimerfiles,probeindices,moleculeindices,numberprobeatoms=GenerateInitialProbeStructure(poltype,missingvdwatomindices)
    obConversion = openbabel.OBConversion()
    for i in range(len(dimerfiles)):
        filenamelist=dimerfiles[i]
        probeindexlist=probeindices[i]
        moleculeindexlist=moleculeindices[i]
        probeatoms=numberprobeatoms[i]
        alloutputfilenames=[]
        prefixarrays=[]
        distancearrays=[]
        checkarray=[]
        for filenameidx in range(len(filenamelist)):
            filename=filenamelist[filenameidx]
            dimermol = openbabel.OBMol()
            probeindex=probeindexlist[filenameidx]
            moleculeindex=moleculeindexlist[filenameidx]
            inFormat = obConversion.FormatFromExt(filename)
            obConversion.SetInFormat(inFormat)
            obConversion.ReadFile(dimermol, filename)
            prefix=filename.replace('.xyz','')
            prefixarrays.append(prefix)
            check=CheckIfFittingCompleted(poltype,prefix)
            checkarray.append(check)
            poltype.comoptfname=prefix+'.com'
            poltype.chkoptfname=prefix+'.chk'
            poltype.fckoptfname=prefix+'.fchk'
            poltype.logoptfname=prefix+'.log'
            poltype.gausoptfname=prefix+'.log'
            optmol = opt.GeometryOptimization(poltype,dimermol,checkbonds=False,modred=False)
            dimeratoms=dimermol.NumAtoms()
            moleculeatoms=dimeratoms-probeatoms
            moleculeatom=optmol.GetAtom(moleculeindex)
            probeatom=optmol.GetAtom(probeindex)
            moleculeatomcoords=np.array([moleculeatom.GetX(),moleculeatom.GetY(),moleculeatom.GetZ()])
            probeatomcoords=np.array([probeatom.GetX(),probeatom.GetY(),probeatom.GetZ()])
            equildistance=np.linalg.norm(probeatomcoords-moleculeatomcoords)
            distarray=np.multiply(equildistance,np.array(array))
            distancearrays.append(distarray)
            outputprefixname=filename.split('.')[0]
            outputxyz=outputprefixname+'_tinker.xyz'
            inputxyz=outputprefixname+'_cartesian.xyz'
            WriteOutCartesianXYZ(poltype,optmol,inputxyz)
            if 'water' in outputxyz:
                waterbool=True
            else:
                waterbool=False
            ConvertProbeDimerXYZToTinkerXYZ(poltype,inputxyz,poltype.xyzoutfile,outputxyz,waterbool,probeatoms)
            filenamearray=MoveDimerAboutMinima(poltype,outputxyz,outputprefixname,moleculeatoms,moleculeindex,probeindex,equildistance,array)
            qmfilenamearray=GenerateSPInputFiles(poltype,filenamearray,poltype.mol,probeatoms)
            outputfilenames=ExecuteSPJobs(poltype,qmfilenamearray,prefix)
            alloutputfilenames.extend(outputfilenames)
        dothefit=False
        for check in checkarray:
            if check==False:
                dothefit=True
        if dothefit==True:
            ReadCounterPoiseAndWriteQMData(poltype,alloutputfilenames)
            vdwtype=poltype.idxtosymclass[moleculeindex]
            vdwtypesarray=[vdwtype]
            initialvdwradius,initialvdwdepth,minvdwradius,maxvdwradius,minvdwdepth,maxvdwdepth=GrabVdwParameters(poltype,vdwtype)
            initialradii=[initialvdwradius]
            initialdepths=[initialvdwdepth]
            minradii=[minvdwradius]
            maxradii=[maxvdwradius]
            mindepths=[minvdwdepth]
            maxdepths=[maxvdwdepth]
            WriteInitialPrmFile(poltype,vdwtypesarray,initialradii,initialdepths,minradii,maxradii,mindepths,maxdepths)
            vdwradius,vdwdepth=VDWOptimizer(poltype)

            for k in range(len(prefixarrays)):
                prefix=prefixarrays[k]
                distarray=distancearrays[k]
                PlotEnergyVsDistance(poltype,distarray,prefix,vdwradius,vdwdepth,vdwtypesarray)
                PlotQMVsMMEnergy(poltype,vdwtypesarray,prefix)
    shutil.copy(poltype.key5fname,'../'+poltype.key5fname)
    os.chdir(poltype.parentdir)
