from __future__ import print_function
from pymol.cgo import *    # get constants
from math import *
from pymol import cmd
import os
import sys
import numpy as np
import re


'''
http://pymolwiki.org/index.php/cgo_arrow

(c) 2013 Thomas Holder, Schrodinger Inc.

License: BSD-2-Clause
'''

from pymol import cmd, cgo, CmdException


def cgo_arrow(atom1='pk1', atom2='pk2', radius=0.5, gap=0.0, hlength=-1, hradius=-1,
              color='blue red', name='',lab=None):
    '''
DESCRIPTION

    Create a CGO arrow between two picked atoms.

ARGUMENTS

    atom1 = string: single atom selection or list of 3 floats {default: pk1}

    atom2 = string: single atom selection or list of 3 floats {default: pk2}

    radius = float: arrow radius {default: 0.5}

    gap = float: gap between arrow tips and the two atoms {default: 0.0}

    hlength = float: length of head

    hradius = float: radius of head

    color = string: one or two color names {default: blue red}

    name = string: name of CGO object
    '''
    from chempy import cpv

    radius, gap = float(radius), float(gap)
    hlength, hradius = float(hlength), float(hradius)

    try:
        color1, color2 = color.split()
    except:
        color1 = color2 = color
    color2=color1
    color1 = list(cmd.get_color_tuple(color1))
    color2 = list(cmd.get_color_tuple(color2))

    def get_coord(v):
        if not isinstance(v, str):
            return v
        if v.startswith('['):
            return cmd.safe_list_eval(v)
        return cmd.get_atom_coords(v)

    xyz1 = get_coord(atom1)
    xyz2 = get_coord(atom2)
    
    normal = cpv.normalize(cpv.sub(xyz1, xyz2))

    

    if gap:
        diff = cpv.scale(normal, gap)
        xyz1 = cpv.sub(xyz1, diff)
        xyz2 = cpv.add(xyz2, diff)
    if not name:
        name = cmd.get_unused_name('cylinder')

    
    name1 = cmd.get_unused_name('arrow')

    
    radius=.05
    x1=xyz1[0]
    y1=xyz1[1]
    z1=xyz1[2]
    x2=xyz2[0]
    y2=xyz2[1]
    z2=xyz2[2]
    r1=color1[0]
    g1=color1[1]
    b1=color1[2]
    r2=color2[0]
    g2=color2[1]
    b2=color2[2]
    
    cmd.load_cgo( [ 9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2 ], name )
    hlength=.5
    xyz3 = cpv.add(-1*np.array(cpv.scale(normal, hlength)), xyz2)
    x3=xyz3[0]
    y3=xyz3[1]
    z3=xyz3[2]
    
    cmd.load_cgo( [ CONE, x2,y2, z2, x3, y3, z3, 2*radius, 0.0, r1, g1, b1, r1, g1, b1, 1.0, 0.0 ], name1 )
    if lab!=None:
        cmd.pseudoatom(pos=[x3,y3,z3],label=lab)


## made by Brandon Walker University of Texas at Austin

def AssumePoltypeDefaults():
    poltypejobpath=ReadPoltypeFilePath()
    os.chdir(poltypejobpath)
    files=os.listdir(os.getcwd())
    for f in files:
        if '.' in f and 'ttt' not in f:
            filesplit=f.split('.')
            suffix=filesplit[1]
            if suffix=='key':
                DefineTinkerKeyFilePath(poltypejobpath+f)
            elif suffix=='key_5':
                DefineTinkerFinalKeyFilePath(poltypejobpath+f)
            elif suffix=='xyz':
                DefineTinkerXYZFilePath(poltypejobpath+f)
            elif suffix=='xyz_2':
                DefineTinkerFinalXYZFilePath(poltypejobpath+f)

            
def DefineTinkerFinalKeyFilePath(keyfilepath):
    tempname='finalkeyfilepath.txt'
    temp=open(tempname,'w')
    temp.write(keyfilepath)
    temp.close()

def DefinePoltypeFilePath(filepath):
    tempname='poltypefilepath.txt'
    temp=open(tempname,'w')
    temp.write(filepath)
    temp.close()


def DefineTinkerKeyFilePath(keyfilepath):
    tempname='keyfilepath.txt'
    temp=open(tempname,'w')
    temp.write(keyfilepath)
    temp.close()



def DefineTinkerFinalXYZFilePath(pathtotinkxyz):
    tempname='pathtotinkfinalxyz.txt'
    temp=open(tempname,'w')
    temp.write(pathtotinkxyz)
    temp.close()

def DefineTinkerXYZFilePath(pathtotinkxyz):
    tempname='pathtotinkxyz.txt'
    temp=open(tempname,'w')
    temp.write(pathtotinkxyz)
    temp.close()
    cmd.load(pathtotinkxyz)
    cmd.show_as('sticks')
    cmd.bg_color("white")
    LabelIndexAndTinkerTypeNumbers()

def ReadTinkerFinalXYZFilePath():
    if os.path.isfile('pathtotinkfinalxyz.txt'):
        temp=open('pathtotinkfinalxyz.txt','r')
        results=temp.readlines()
        temp.close()
        firstline=results[0]
        linesplit=firstline.split()
        return linesplit[0]
    else:
        print('Tinker Final XYZ file path has not been defined')

def ReadTinkerXYZFilePath():
    if os.path.isfile('pathtotinkxyz.txt'):
        temp=open('pathtotinkxyz.txt','r')
        results=temp.readlines()
        temp.close()
        firstline=results[0]
        linesplit=firstline.split()
        return linesplit[0]
    else:
        print('Tinker XYZ file path has not been defined')

def ReadTinkerFinalKeyFilePath():
    if os.path.isfile('finalkeyfilepath.txt'):
        temp=open('finalkeyfilepath.txt','r')
        results=temp.readlines()
        temp.close()
        firstline=results[0]
        linesplit=firstline.split()
        return linesplit[0]
    else:
        print('Final Key file path has not been defined')

def ReadPoltypeFilePath():
    if os.path.isfile('poltypefilepath.txt'):
        temp=open('poltypefilepath.txt','r')
        results=temp.readlines()
        temp.close()
        firstline=results[0]
        linesplit=firstline.split()
        return linesplit[0]
    else:
        print('poltype file path has not been defined')

def ReadTinkerKeyFilePath():
    if os.path.isfile('keyfilepath.txt'):
        temp=open('keyfilepath.txt','r')
        results=temp.readlines()
        temp.close()
        firstline=results[0]
        linesplit=firstline.split()
        return linesplit[0]
    else:
        print('Key file path has not been defined')


def CallSubprocess(cmdstr):
    os.system(cmdstr)


            

def DefineTinkerBINPath(binpath):
    tempname='binpath.txt'
    temp=open(tempname,'w')
    temp.write(binpath)
    temp.close()

def ReadTinkerBINPath():
    if os.path.isfile('binpath.txt'):
        temp=open('binpath.txt','r')
        results=temp.readlines()
        temp.close()
        firstline=results[0]
        linesplit=firstline.split()
        return linesplit[0]
    else:
        print('binary path has not been defined')

def GrabTinkerExecutablePath(string):
    curpath=os.getcwd()
    binpath=ReadTinkerBINPath()
    os.chdir(binpath)
    files=os.listdir(os.getcwd())
    found=False
    for f in files:
        if string in f:
            executable=f
            found=True
    if found==True:
        return binpath+'/'+executable
    else:
        print(executablestring+' not found in '+binpath)

def PrintTinkerTotalCharge():
    alzpath=GrabTinkerExecutablePath('analyze')
    keyfilepath=ReadTinkerFinalKeyFilePath()
    xyzfilepath=ReadTinkerFinalXYZFilePath()
    alzoutfilepath=os.getcwd()+'/'+'alz.out'
    alzcmd=alzpath+' '+xyzfilepath+' '+'-k'+' '+keyfilepath+' '+'m'+' '+'>'+' '+alzoutfilepath
    CallSubprocess(alzcmd)
    temp=open(alzoutfilepath,'r')
    results=temp.readlines()
    for line in results:
        if 'Total Electric Charge :' in line:
            print(line)


def PrintTinkerComponentEnergy():
    alzpath=GrabTinkerExecutablePath('analyze')
    keyfilepath=ReadTinkerFinalKeyFilePath()
    xyzfilepath=ReadTinkerFinalXYZFilePath()
    alzoutfilepath=os.getcwd()+'/'+'alz.out'
    alzcmd=alzpath+' '+xyzfilepath+' '+'-k'+' '+keyfilepath+' '+'e'+' '+'>'+' '+alzoutfilepath
    CallSubprocess(alzcmd)
    temp=open(alzoutfilepath,'r')
    results=temp.readlines()
    found=False
    for line in results:
        if 'Total Potential Energy :' in line or found==True:
            print(line)
            found=True
       

def GenerateIndexToCoordsDic():
    pathtotinkxyz=ReadTinkerFinalXYZFilePath()
    temp=open(pathtotinkxyz,'r')
    results=temp.readlines()
    temp.close()
    indextocoords={}
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if lineidx!=0 and len(linesplit)>=2:
            index=linesplit[0]
            x=float(linesplit[2])
            y=float(linesplit[3])
            z=float(linesplit[4])
            vec=[x,y,z]
            indextocoords[index]=vec
    return indextocoords
            
def GenerateIndexToFrameUnitVectors(indextocoords,indextotypenum,indextoconn,indextoframedic):
    #print(indextoframedic) 
    indextoframeunitvectors={}
    indextoframeunitvectors['x']={}
    indextoframeunitvectors['y']={}
    indextoframeunitvectors['z']={}
    for index in indextotypenum.keys():
        typenum=indextotypenum[index]
        coords=np.array(indextocoords[index])
        framedef=indextoframedic['framedef'][index]
        neighbs=indextoframedic['neighbors'][index]
        if framedef=='undefined':
            indextoframeunitvectors['x'][index]=[0,0,0]
            indextoframeunitvectors['y'][index]=[0,0,0]
            indextoframeunitvectors['z'][index]=[0,0,0]
        elif framedef=='z-only':

            mzindex=GrabFirstNeighborIndex(index)

            mzcoords=np.array(indextocoords[mzindex])
            z=np.array(mzcoords-coords)
            mag=np.linalg.norm(z)
            uz=z/mag
            indextoframeunitvectors['x'][index]=[0,0,0]
            indextoframeunitvectors['y'][index]=[0,0,0]
            indextoframeunitvectors['z'][index]=uz
        elif framedef=='z-then-x':
            #print('index ',index,'this is z then x') 
            mzindex=GrabFirstNeighborIndex(index)
            mxindex=GrabSecondNeighborIndex(index)
            mzcoords=np.array(indextocoords[mzindex])
            z=np.array(mzcoords-coords)
            mag=np.linalg.norm(z)
            uz=z/mag
            mxcoords=np.array(indextocoords[mxindex])
            v=np.array(mxcoords-coords)
            uv=v/np.linalg.norm(v)
            uy=np.cross(uv,uz)
            ux=np.cross(uz,uy)
            indextoframeunitvectors['x'][index]=ux
            indextoframeunitvectors['y'][index]=uy
            indextoframeunitvectors['z'][index]=uz

        elif framedef=='bisector':
            mzindex=GrabFirstNeighborIndex(index)
            mxindex=GrabSecondNeighborIndex(index)
            mzcoords=np.array(indextocoords[mzindex])
            mxcoords=np.array(indextocoords[mxindex])
            v1=np.array(mxcoords-coords)
            v2=np.array(mzcoords-coords)
            v3=v2-v1
            disp=v3*.5
            z=disp+v1
            mag=np.linalg.norm(z)
            uz=z/mag
            indextoframeunitvectors['x'][index]=[0,0,0]
            indextoframeunitvectors['y'][index]=[0,0,0]
            indextoframeunitvectors['z'][index]=uz


        elif framedef=='z-then-bisector':

            mzindex=GrabFirstNeighborIndex(index)
            mxindex=GrabSecondNeighborIndex(index)
            m1index=GrabThirdNeighborIndex(index)

            mzcoords=np.array(indextocoords[mzindex])
            z=np.array(mzcoords-coords)
            mag=np.linalg.norm(z)
            uz=z/mag
            mxcoords=np.array(indextocoords[mxindex])
            m1coords=np.array(indextocoords[m1index])
            v1=np.array(mxcoords-coords)
            v2=np.array(m1coords-coords)
            v3=v2-v1
            disp=v3*.5
            xtemp=disp+v1
            magx=np.linalg.norm(xtemp)
            uxtemp=xtemp/magx
            y=np.cross(uz,uxtemp)
            uy=y/np.linalg.norm(y)
            ux=np.cross(uy,uz)
            indextoframeunitvectors['x'][index]=ux
            indextoframeunitvectors['y'][index]=uy
            indextoframeunitvectors['z'][index]=uz

        elif framedef=='trisector':
            mzindex=GrabFirstNeighborIndex(index)
            mxindex=GrabSecondNeighborIndex(index)
            m1index=GrabThirdNeighborIndex(index)
            mzcoords=indextocoords[mzindex]
            mxcoords=indextocoords[mxindex]
            mycoords=indextocoords[m1index]
            v1=np.array(mxcoords-coords)
            v2=np.array(mzcoords-coords)
            v3=np.array(mycoords-coords)
            diff1=v1-v2
            udiff1=diff1/np.linalg.norm(diff1)
            diff2=v2-v3
            udiff2=diff2/np.linalg.norm(diff2)
            z=np.cross(udiff1,udiff2) 
            indextoframeunitvectors['x'][index]=[0,0,0]
            indextoframeunitvectors['y'][index]=[0,0,0]
            indextoframeunitvectors['z'][index]=z

    return indextoframeunitvectors            

def GenerateIndexToNewXYZStructures(indextoframeunitvectors,indextocoords):
    indextozstructurecoords={}
    indextoystructurecoords={}
    indextoxstructurecoords={}
    for index in indextocoords.keys():
        ux=indextoframeunitvectors['x'][index]
        uy=indextoframeunitvectors['y'][index]
        uz=indextoframeunitvectors['z'][index]
        coords=indextocoords[index]
        xcoords=coords+ux
        ycoords=coords+uy
        zcoords=coords+uz
        indextoxstructurecoords[index]=xcoords
        indextoystructurecoords[index]=ycoords
        indextozstructurecoords[index]=zcoords
    return indextoxstructurecoords,indextoystructurecoords,indextozstructurecoords


def GrabFirstNeighborIndex(index):
    keyfilepath=ReadTinkerKeyFilePath()
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'multipole' in line:
            linesplit=line.split()
            atomidx=linesplit[1]
            if atomidx==index:
                firstneighb=linesplit[2].replace('-','')
                return firstneighb
    
     
def GrabSecondNeighborIndex(index):
    keyfilepath=ReadTinkerKeyFilePath()
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'multipole' in line:
            linesplit=line.split()
            atomidx=linesplit[1]
            if atomidx==index:
                secondneighb=linesplit[3].replace('-','')
                return secondneighb
    

def GrabThirdNeighborIndex(index):
    keyfilepath=ReadTinkerKeyFilePath()
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'multipole' in line:
            linesplit=line.split()
            atomidx=linesplit[1]
            if atomidx==index:
                thirdneighb=linesplit[4].replace('-','')
                return thirdneighb

def ShowTinkerMultipoleFrame(index):
    pathtotinkxyz=ReadTinkerFinalXYZFilePath()
    indextoconn=GenerateIndexToConnectivityDic()
    indextocoords=GenerateIndexToCoordsDic()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    indextoframedic=GenerateIndexToFrameDic() # try this first else, have to do it by type numbers, if using poledit, this should take priority over typenumbers
    indextoframeunitvectors=GenerateIndexToFrameUnitVectors(indextocoords,indextotypenum,indextoconn,indextoframedic)
    indextoxstructurecoords,indextoystructurecoords,indextozstructurecoords=GenerateIndexToNewXYZStructures(indextoframeunitvectors,indextocoords) 
    xcoords=indextoxstructurecoords[index]
    ycoords=indextoystructurecoords[index]
    zcoords=np.array(indextozstructurecoords[index])
    framedef=indextoframedic['framedef'][index]
    print('Frame definition: '+framedef)
    PrintTinkerMultipoleParameters(index)
    coords=indextocoords[index]
    cgo_arrow(coords,xcoords,color='red',lab='x')
    cgo_arrow(coords,ycoords,color='blue',lab='y')
    cgo_arrow(coords,zcoords,color='green',lab='z')
    return xcoords,ycoords,zcoords,coords

    

            
def GenerateIndexToConnectivityDic():
    indextoconn={}
    pathtotinkxyz=ReadTinkerFinalXYZFilePath()
    temp=open(pathtotinkxyz,'r')
    results=temp.readlines()
    temp.close()
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if lineidx!=0 and len(linesplit)>=2:
            index=linesplit[0]
            neighbs=linesplit[6:]
            indextoconn[index]=neighbs
    return indextoconn

def GenerateIndexToTinkerTypeDic():
    pathtotinkxyz=ReadTinkerFinalXYZFilePath()
    temp=open(pathtotinkxyz,'r')
    results=temp.readlines()
    temp.close()
    indextotypenum={}
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if lineidx!=0 and len(linesplit)>=2:
            index=linesplit[0]
            type=linesplit[5]
            indextotypenum[index]=type
    return indextotypenum


def GenerateIndexToFrameDic():
    indextoframedic={}
    keyfilepath=ReadTinkerKeyFilePath()
    prmfilepath=None
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    indextoframedic['framedef']={}
    indextoframedic['neighbors']={}
    indextotypenum=GenerateIndexToTinkerTypeDic()
    for line in results:
        if 'multipole' in line and '#' not in line:
            linesplit=line.split()
            index=linesplit[1].replace('-','')
            if index not in indextotypenum.values(): # if not a type number
                if ((len(linesplit)==3 and '0' not in linesplit) or (len(linesplit)==5 and '0'==linesplit[2] and '0'==linesplit[3])):
                    framedef='undefined'
                    neighbs=[]
                if ((len(linesplit)==4 and '0' not in linesplit) or (len(linesplit)==5 and '0' in linesplit)) and '-' not in linesplit[2:-1]:
                    framedef='z-only'
                    neighbs=[linesplit[2]]
                if (len(linesplit)==5 and '0' not in linesplit) and '-' not in ' '.join(linesplit[2:4]):
                    framedef='z-then-x'
                    neighbs=[linesplit[2].replace('-',''),linesplit[3].replace('-','')]
                if (len(linesplit)==5 and '-' in linesplit[2] and '-' not in linesplit[3]) or (len(linesplit)==5 and '-' in linesplit[3] and '-' not in linesplit[2]) or (len(linesplit)==5 and '-' in linesplit[2] and '-' in linesplit[3]):
                    framedef='bisector'
                    neighbs=[linesplit[2].replace('-',''),linesplit[3].replace('-','')]
                if (len(linesplit)==6 and '-' in linesplit[3] and '-' in linesplit[4] and '-' not in linesplit[2]):
                    framedef='z-then-bisector'
                    neighbs=[linesplit[2],linesplit[3].replace('-',''),linesplit[4].replace('-','')]
                if (len(linesplit)==6 and '-' in linesplit[2] and '-' in linesplit[3] and '-' in linesplit[4]):
                    framedef='trisector'
                    neighbs=[linesplit[2].replace('-',''),linesplit[3].replace('-',''),linesplit[4].replace('-','')]
                
                indextoframedic['framedef'][index]=framedef
                indextoframedic['neighbors'][index]=neighbs
    


    return indextoframedic
                

def GrabParameterFilePathFromTinkerKey():
    keyfilepath=ReadTinkerFinalKeyFilePath()
    prmfilepath=None
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'parameters' in line:
            linesplit=line.split()
            prmfilepath=linesplit[1]
            return prmfilepath
    if prmfilepath==None:
        print('Parameter file path not found in keyfile '+keyfilepath)

def GenerateTinkerTypeToTinkerClassDic():
    keyfilepath=ReadTinkerFinalKeyFilePath()
    pathtotinkxyz=ReadTinkerFinalXYZFilePath()
    atomtypetoatomclass={}
    pathtotinkxyz=ReadTinkerFinalXYZFilePath()
    indextotypedic=GenerateIndexToTinkerTypeDic()
    prmfilepath=GrabParameterFilePathFromTinkerKey()
    temp=open(prmfilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'atom' in line:
            linesplit=line.split()
            typenum=linesplit[1]
            classnum=linesplit[2]
            if typenum in indextotypedic.values():
                atomtypetoatomclass[typenum]=classnum
  
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'atom' in line:
            linesplit=line.split()
            typenum=linesplit[1]
            classnum=linesplit[2]
            if typenum in indextotypedic.values():
                atomtypetoatomclass[typenum]=classnum
    return atomtypetoatomclass
                


def LabelIndexAndTinkerTypeNumbers():
    indextotypenum=GenerateIndexToTinkerTypeDic()
    for index in indextotypenum.keys():
        typenum=indextotypenum[index]
        lab=index+','+typenum
        cmd.select("index "+index)
        cmd.label("sele",lab)

def LabelTinkerTypeNumbers():
    indextotypenum=GenerateIndexToTinkerTypeDic()
    for index in indextotypenum.keys():
        typenum=indextotypenum[index]
        lab=typenum
        cmd.select("index "+index)
        cmd.label("sele",lab)

def LabelIndexAndTinkerTypeAndClassNumbers():
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    for index in indextotypenum.keys():
        typenum=indextotypenum[index]
        classnum=atomtypetoatomclass[typenum]
        lab=index+','+typenum+','+classnum
        cmd.select("index "+index)
        cmd.label("sele",lab)

def LabelTinkerClassNumbers():
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    for index in indextotypenum.keys():
        typenum=indextotypenum[index]
        classnum=atomtypetoatomclass[typenum]
        lab=classnum
        cmd.select("index "+index)
        cmd.label("sele",lab)

def GrabPolarizationParameters(typenumber,prmfilepath):
    polline=None
    temp=open(prmfilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'polarize' in line:
            linesplit=line.split()
            if len(linesplit)>=5:
                poltype=linesplit[1]
                if poltype==typenumber:
                    polline=line
                    return polline
    return polline

def GrabMultipoleParameters(typenumber,prmfilepath,atomindex=None): # if index is defined, check for this frame as priority over typenumber( if poledit was used)
    indextotypenum=GenerateIndexToTinkerTypeDic()
    print('indextotypenum',indextotypenum)
    indextocoords=GenerateIndexToCoordsDic()
    mpllines=[]
    temp=open(prmfilepath,'r')
    results=temp.readlines()
    temp.close()
    for lineidx in range(len(results)):
        mpllines=CheckMultipoleResult(results,typenumber,lineidx,mpllines,indextotypenum)

    return mpllines

 
def GrabAllIndexesForType(typenum,indextotypenum):
    print(' typenum ',typenum)
    ls=[]
    for index in indextotypenum.keys():
        typeidx=indextotypenum[index]
        print('typeidx ',typeidx)
        if typeidx==typenum:
            print(' typeidx ',typeidx,'type number ',typenumber,' index ',index)
            ls.append(index)
    return ls


 
def CheckMultipoleResult(results,typenumber,lineidx,mpllines,indextotypenum):
    line=results[lineidx]
    if 'multipole' in line:
        linesplit=line.split()
        if len(linesplit)>=5:
            poltype=linesplit[1]
            if poltype==typenumber:
                    
                currentneighbs=linesplit[2:-1]
                foundall=True

                for neighb in currentneighbs:
                    newneighb=neighb.replace('-','')
                    
                    if newneighb not in indextotypenum.values():
                                                   
                        foundall=False
                if foundall==True:
                    if results[lineidx] not in mpllines:
                        mpllines.append(results[lineidx])
                        mpllines.append(results[lineidx+1])
                        mpllines.append(results[lineidx+2])
                        mpllines.append(results[lineidx+3])
                        mpllines.append(results[lineidx+4])
                    return mpllines
    return mpllines

   

def PrintTinkerMultipoleParameters(atomindex):
    pathtotinkxyz=ReadTinkerFinalXYZFilePath()
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    typenumber=indextotypenum[atomindex]
    prmfilepath=GrabParameterFilePathFromTinkerKey()
    indextoconn=GenerateIndexToConnectivityDic()
    neighbs=[]
    conn=indextoconn[atomindex]
    print(keyfilepath)
    mpllines=GrabMultipoleParameters(typenumber,keyfilepath,atomindex)
    if len(mpllines)==0:
         mpllines=GrabMultipoleParameters(typenumber,keyfilepath,atomindex=None)
    if len(mpllines)==0:
         mpllines=GrabMultipoleParameters(typenumber,prmfilepath,atomindex=None)
    if len(mpllines)==0:
        print('Parameters Not Found')
    else:
        for line in mpllines:
            print(line)
    return mpllines

def GrabVdwParameters(classnumber,prmfilepath):
    vdwline=None
    temp=open(prmfilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if 'vdw' in line and 'vdwpr' not in line and len(linesplit)>=4:
            cls=linesplit[1]
            if cls==classnumber:
                vdwline=line
                return vdwline
    return vdwline

def GenerateNewKeyName():
    keyfilepath=ReadTinkerFinalKeyFilePath()
    curpath=os.getcwd()
    dir_path = os.path.dirname(os.path.realpath(keyfilepath))
    os.chdir(dir_path)
    files=os.listdir(os.getcwd())
    lastnumber=1
    for f in files:
        if 'key_' in f:
            idx=f.find('key_')
            numstartidx=idx+4
            num=int(f[numstartidx:])
            if num>lastnumber:
                lastnumber=num
    if 'key_' in keyfilepath:
        idx=keyfilepath.find('key_')
        finalidx=idx+2
        string=keyfilepath[:finalidx+1]
    else:
        string=keyfilepath
    
    newkeyname=string.replace('.key','.key_'+str(lastnumber+1))
    print('New keyname is '+newkeyname)
    return newkeyname

def ChangeTinkerVdwParameters(atomindex,prm1,prm2):
    indextotypenum=GenerateIndexToTinkerTypeDic()
    typenumber=indextotypenum[atomindex]
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    classnumber=atomtypetoatomclass[typenumber]    
    keyfilepath=ReadTinkerFinalKeyFilePath()
    newkeyname=GenerateNewKeyName()
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    temp=open(newkeyname,'w')
    foundprm=False
    for line in results:
        linesplit=line.split()
        if 'vdw' in line and 'vdwpr' not in line and len(linesplit)>=4:
            cls=linesplit[1]
            if cls==classnumber:
                foundprm=True
                oldprm1=linesplit[2]
                oldprm2=linesplit[3]
                linesplit[2]=prm1
                linesplit[3]=prm2
                newline=' '.join(linesplit)
                temp.write(newline+'\n')
            else:
               temp.write(line)
        else:
            temp.write(line)
    temp.close()
    if foundprm==False:
        print('Parameters are not in keyfile, need to call AddTinkerVdwParameters to add to keyfile')

def AddTinkerVdwParameters(atomindex,prm1,prm2):
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    typenumber=indextotypenum[atomindex]
    classnumber=atomtypetoatomclass[typenumber]
    vdwline='vdw '+classnumber+' '+prm1+' '+prm2+' '+'\n'
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    newkeyname=GenerateNewKeyName()
    temp=open(newkeyname,'w')
    foundfirst=False
    wrote=False
    for lineidx in range(len(results)):
        line=results[lineidx]
        if 'vdw' in line and 'lambda' not in line and 'cutoff' not in line:
            foundfirst=True
            temp.write(line)
            temp.flush()
            os.fsync(temp.fileno())
            sys.stdout.flush()
        else:
            if foundfirst==True and wrote==False:
                temp.write(vdwline)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
                wrote=True
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            else:
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
    temp.close()

def AddTinkerPolarizationParameters(atomindex,prm1,prm2,neighboridxs):
    neighbidxlist=neighboridxs.split()
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    typenumber=indextotypenum[atomindex]
    neighbtypes=[]
    for idx in neighbidxlist:
        typenum=indextotypenum[idx]
        neighbtypes.append(typenum)
    polline='polarize '+typenumber+' '+prm1+' '+prm2+' '
    for neighbtype in neighbtypes:
        polline+=neighbtype+' '
    polline+='\n'
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    newkeyname=GenerateNewKeyName()
    temp=open(newkeyname,'w')
    foundfirst=False
    wrote=False
    for lineidx in range(len(results)):
        line=results[lineidx]
        if 'polarize' in line:
            foundfirst=True
            temp.write(line)
            temp.flush()
            os.fsync(temp.fileno())
            sys.stdout.flush()
        else:
            if foundfirst==True and wrote==False:
                temp.write(polline)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
                wrote=True
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            else:
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
    temp.close()
    
def ChangeTinkerPolarizationParameters(atomindex,prm1,prm2):
    indextotypenum=GenerateIndexToTinkerTypeDic()
    typenumber=indextotypenum[atomindex]
    keyfilepath=ReadTinkerFinalKeyFilePath()
    newkeyname=GenerateNewKeyName()
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    temp=open(newkeyname,'w')
    foundprm=False
    for line in results:
        linesplit=line.split()
        if 'polarize' in line and len(linesplit)>=5:
            typenum=linesplit[1]
            print(typenum,typenumber)
            if typenum==typenumber:
                foundprm=True
                oldprm1=linesplit[2]
                oldprm2=linesplit[3]
                linesplit[2]=prm1
                linesplit[3]=prm2
                newline=' '.join(linesplit)
                temp.write(newline+'\n')
            else:
               temp.write(line)
        else:
            temp.write(line)
    temp.close()
    if foundprm==False:
        print('Parameters are not in keyfile, need to call AddTinkerPolarizationParameters to add to keyfile')

def ChangeTinkerBondParameters(atomindex1,atomindex2,prm1,prm2):
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    typenumber1=indextotypenum[atomindex1]
    typenumber2=indextotypenum[atomindex2]
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    classnumber1=atomtypetoatomclass[typenumber1]
    classnumber2=atomtypetoatomclass[typenumber2]
    newkeyname=GenerateNewKeyName()
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    temp=open(newkeyname,'w')
    foundprm=False
    for line in results:
        linesplit=line.split()
        if 'bond' in line and len(linesplit)>=5:
            if classnumber1 in linesplit and classnumber2 in linesplit:
                foundprm=True
                oldprm1=linesplit[3]
                oldprm2=linesplit[4]
                linesplit[3]=prm1
                linesplit[4]=prm2
                newline=' '.join(linesplit)
                temp.write(newline+'\n')
            else:
               temp.write(line)
        else:
            temp.write(line)
    temp.close()
    if foundprm==False:
        print('Parameters are not in keyfile, need to call AddTinkerBondParameters to add to keyfile')

def ChangeTinkerAngleParameters(atomindex1,atomindex2,atomindex3,prm1,prm2):
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    typenumber1=indextotypenum[atomindex1]
    typenumber2=indextotypenum[atomindex2]
    typenumber3=indextotypenum[atomindex3]
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    classnumber1=atomtypetoatomclass[typenumber1]
    classnumber2=atomtypetoatomclass[typenumber2]
    classnumber3=atomtypetoatomclass[typenumber3]
    newkeyname=GenerateNewKeyName()
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    temp=open(newkeyname,'w')
    foundprm=False
    for line in results:
        linesplit=line.split()
        if 'angle' in line and len(linesplit)>=6:
            if classnumber1 in linesplit and classnumber2 in linesplit and classnumber3 in linesplit:
                foundprm=True
                oldprm1=linesplit[4]
                oldprm2=linesplit[5]
                linesplit[4]=prm1
                linesplit[5]=prm2
                newline=' '.join(linesplit)
                temp.write(newline+'\n')
            else:
               temp.write(line)
        else:
            temp.write(line)
    temp.close()
    if foundprm==False:
        print('Parameters are not in keyfile, need to call AddTinkerAngleParameters to add to keyfile')

def ChangeTinkerTorsionParameters(atomindex1,atomindex2,atomindex3,atomindex4,prm1,prm2,prm3):
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    typenumber1=indextotypenum[atomindex1]
    typenumber2=indextotypenum[atomindex2]
    typenumber3=indextotypenum[atomindex3]
    typenumber4=indextotypenum[atomindex4]
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    classnumber1=atomtypetoatomclass[typenumber1]
    classnumber2=atomtypetoatomclass[typenumber2]
    classnumber3=atomtypetoatomclass[typenumber3]
    classnumber4=atomtypetoatomclass[typenumber4]
    newkeyname=GenerateNewKeyName()
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    temp=open(newkeyname,'w')
    foundprm=False
    for line in results:
        linesplit=line.split()
        if 'torsion' in line and len(linesplit)>=6:
            if classnumber1 in linesplit and classnumber2 in linesplit and classnumber3 in linesplit and classnumber4 in linesplit:
                foundprm=True
                oldprm1=linesplit[5]
                oldprm2=linesplit[8]
                oldprm3=linesplit[11]
                linesplit[5]=prm1
                linesplit[8]=prm2
                linesplit[11]=prm3
                newline=' '.join(linesplit)
                temp.write(newline+'\n')
            else:
               temp.write(line)
        else:
            temp.write(line)
    temp.close()
    if foundprm==False:
        print('Parameters are not in keyfile, need to call AddTinkerTorsionParameters to add to keyfile')

def ChangeTinkerStretchBendParameters(atomindex1,atomindex2,atomindex3,prm1,prm2):
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    typenumber1=indextotypenum[atomindex1]
    typenumber2=indextotypenum[atomindex2]
    typenumber3=indextotypenum[atomindex3]
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    classnumber1=atomtypetoatomclass[typenumber1]
    classnumber2=atomtypetoatomclass[typenumber2]
    classnumber3=atomtypetoatomclass[typenumber3]
    newkeyname=GenerateNewKeyName()
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    temp=open(newkeyname,'w')
    foundprm=False
    for line in results:
        linesplit=line.split()
        if 'strbnd' in line and len(linesplit)>=6:
            if classnumber1 in linesplit and classnumber2 in linesplit and classnumber3 in linesplit:
                foundprm=True
                oldprm1=linesplit[4]
                oldprm2=linesplit[5]
                linesplit[4]=prm1
                linesplit[5]=prm2
                newline=' '.join(linesplit)
                temp.write(newline+'\n')
            else:
               temp.write(line)
        else:
            temp.write(line)
    temp.close()
    if foundprm==False:
        print('Parameters are not in keyfile, need to call AddTinkerStretchBendParameters to add to keyfile')

def AddTinkerBondParameters(atomindex1,atomindex2,prm1,prm2):
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    typenumber1=indextotypenum[atomindex1]
    typenumber2=indextotypenum[atomindex2]
    classnumber1=atomtypetoatomclass[typenumber1]
    classnumber2=atomtypetoatomclass[typenumber2]
    bondline='bond '+classnumber1+' '+classnumber2+' '+prm1+' '+prm2+'\n'
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    newkeyname=GenerateNewKeyName()
    temp=open(newkeyname,'w')
    foundfirst=False
    wrote=False
    for lineidx in range(len(results)):
        line=results[lineidx]
        if 'bond' in line:
            foundfirst=True
            temp.write(line)
            temp.flush()
            os.fsync(temp.fileno())
            sys.stdout.flush()
        else:
            if foundfirst==True and wrote==False:
                temp.write(bondline)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
                wrote=True
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            else:
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
    temp.close()

def AddTinkerAngleParameters(atomindex1,atomindex2,atomindex3,prm1,prm2):
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    typenumber1=indextotypenum[atomindex1]
    typenumber2=indextotypenum[atomindex2]
    typenumber3=indextotypenum[atomindex3]
    classnumber1=atomtypetoatomclass[typenumber1]
    classnumber2=atomtypetoatomclass[typenumber2]
    classnumber3=atomtypetoatomclass[typenumber3]
    angleline='angle '+classnumber1+' '+classnumber2+' '+classnumber3+' '+prm1+' '+prm2+'\n'
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    newkeyname=GenerateNewKeyName()
    temp=open(newkeyname,'w')
    foundfirst=False
    wrote=False
    for lineidx in range(len(results)):
        line=results[lineidx]
        if 'angle' in line:
            foundfirst=True
            temp.write(line)
            temp.flush()
            os.fsync(temp.fileno())
            sys.stdout.flush()
        else:
            if foundfirst==True and wrote==False:
                temp.write(angleline)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
                wrote=True
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            else:
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
    temp.close()

def AddTinkerStretchBendParameters(atomindex1,atomindex2,atomindex3,prm1,prm2):
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    typenumber1=indextotypenum[atomindex1]
    typenumber2=indextotypenum[atomindex2]
    typenumber3=indextotypenum[atomindex3]
    classnumber1=atomtypetoatomclass[typenumber1]
    classnumber2=atomtypetoatomclass[typenumber2]
    classnumber3=atomtypetoatomclass[typenumber3]
    strbndline='strbnd '+classnumber1+' '+classnumber2+' '+classnumber3+' '+prm1+' '+prm2+'\n'
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    newkeyname=GenerateNewKeyName()
    temp=open(newkeyname,'w')
    foundfirst=False
    wrote=False
    for lineidx in range(len(results)):
        line=results[lineidx]
        if 'strbnd' in line:
            foundfirst=True
            temp.write(line)
            temp.flush()
            os.fsync(temp.fileno())
            sys.stdout.flush()
        else:
            if foundfirst==True and wrote==False:
                temp.write(strbndline)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
                wrote=True
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            else:
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
    temp.close()

def AddTinkerTorsionParameters(atomindex1,atomindex2,atomindex3,atomindex4,prm1,prm2,prm3):
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    typenumber1=indextotypenum[atomindex1]
    typenumber2=indextotypenum[atomindex2]
    typenumber3=indextotypenum[atomindex3]
    typenumber4=indextotypenum[atomindex4]
    classnumber1=atomtypetoatomclass[typenumber1]
    classnumber2=atomtypetoatomclass[typenumber2]
    classnumber3=atomtypetoatomclass[typenumber3]
    classnumber4=atomtypetoatomclass[typenumber4]
    torsionline='torsion '+classnumber1+' '+classnumber2+' '+classnumber3+' '+classnumber4+' '+prm1+' '+'0.0 1'+' '+prm2+' '+'180.0 2'+' '+prm3+' '+'0.0 3'+'\n'
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    newkeyname=GenerateNewKeyName()
    temp=open(newkeyname,'w')
    foundfirst=False
    wrote=False
    for lineidx in range(len(results)):
        line=results[lineidx]
        if 'torsion' in line:
            foundfirst=True
            temp.write(line)
            temp.flush()
            os.fsync(temp.fileno())
            sys.stdout.flush()
        else:
            if foundfirst==True and wrote==False:
                temp.write(torsionline)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
                wrote=True
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
            else:
                temp.write(line)
                temp.flush()
                os.fsync(temp.fileno())
                sys.stdout.flush()
    temp.close()


def GrabDipoleVector(atomindex):
    indextotypenum=GenerateIndexToTinkerTypeDic()
    typenum=indextotypenum[atomindex]
    keyfilepath=ReadTinkerFinalKeyFilePath()
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    dipolevec=[]
    for lineidx in range(len(results)):
        line=results[lineidx]
        if 'multipole' in line:
            linesplit=line.split()
            atomtype=linesplit[1]
            if atomtype==str(typenum):
                dipoleline=results[lineidx+1]
                dipolelinesplit=dipoleline.split()
                dipolevec.append(float(dipolelinesplit[0]))
                dipolevec.append(float(dipolelinesplit[1]))
                dipolevec.append(float(dipolelinesplit[2]))
    return dipolevec

def IsQxzComponentZero(atomindex):
    iszero=True
    indextotypenum=GenerateIndexToTinkerTypeDic()
    typenum=indextotypenum[atomindex]
    keyfilepath=ReadTinkerFinalKeyFilePath()
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    dipolevec=[]
    for lineidx in range(len(results)):
        line=results[lineidx]
        if 'multipole' in line:
            linesplit=line.split()
            atomtype=linesplit[1]
            if atomtype==str(typenum):
                lineneeded=results[lineidx+4]
                lineneededsplit=lineneeded.split()
                if float(lineneededsplit[0])!=0:
                    iszero=False
    return iszero


def ShowTinkerQuadrupoleCharges(atomindex):
    iszero=IsQxzComponentZero(atomindex)
    xcoords,ycoords,zcoords,coords=ShowTinkerMultipoleFrame(atomindex) 
    r1,r2,r1name,r2name,uvec=ShowTinkerDipoleCharges(atomindex)
    cmd.delete(r1name)
    cmd.delete(r2name)
    r1disp=r1-coords
    monopole=GrabMonopoleCharge(atomindex)
    if iszero==False:
        disp=np.cross(ycoords,uvec)
        newr1=r1+.5*disp # negative charge
        newr2=r1-.5*disp # positive charge
        r3=newr2-2*r1disp # negative charge
        r4=newr1-2*r1disp # positive charge
        newr1name=str(atomindex)+'_'+'FirstNegativeQuadrupoleCharge'
        newr2name=str(atomindex)+'_'+'FirstPositiveQuadrupoleCharge'
        r3name=str(atomindex)+'_'+'SecondNegativeQuadrupoleCharge'
        r4name=str(atomindex)+'_'+'SecondPositiveQuadrupoleCharge'
        newr1label='-'+str(monopole)
        newr2label=str(monopole)
        r3label='-'+str(monopole)
        r4label=str(monopole)
        newr1color='blue'
        newr2color='red'
        r3color='blue'
        r4color='red'
        radius=.05
        cmd.pseudoatom(newr1name,name=newr1name,b=radius,color=newr1color,label=newr1label,pos=list(newr1))
        cmd.pseudoatom(newr2name,name=newr2name,b=radius,color=newr2color,label=newr2label,pos=list(newr2))
        cmd.pseudoatom(r3name,name=r3name,b=radius,color=r3color,label=r3label,pos=list(r3))
        cmd.pseudoatom(r4name,name=r4name,b=radius,color=r4color,label=r4label,pos=list(r4))
        r1string='nb_spheres, '+'('+'name'+' '+newr1name+')'
        r2string='nb_spheres, '+'('+'name'+' '+newr2name+')'
        r3string='nb_spheres, '+'('+'name'+' '+r3name+')'
        r4string='nb_spheres, '+'('+'name'+' '+r4name+')'
        cmd.show(r1string)
        cmd.show(r2string)
        cmd.show(r3string)
        cmd.show(r4string)
    else: # then axial quadrupole
        newr1=coords+2*r1disp # negative q
        newr2=coords-2*r1disp # negative q
        r3=coords # positive 2q
        newr1name=str(atomindex)+'_'+'FirstNegativeQuadrupoleCharge'
        newr2name=str(atomindex)+'_'+'SecondNegativeQuadrupoleCharge'
        r3name=str(atomindex)+'_'+'PositiveQuadrupoleCharge'
        newr1label='-'+str(monopole)
        newr2label='-'+str(monopole)
        r3label=str(2*monopole)
        newr1color='blue'
        newr2color='blue'
        r3color='red'
        radius=.05
        cmd.pseudoatom(newr1name,name=newr1name,b=radius,color=newr1color,label=newr1label,pos=list(newr1))
        cmd.pseudoatom(newr2name,name=newr2name,b=radius,color=newr2color,label=newr2label,pos=list(newr2))
        cmd.pseudoatom(r3name,name=r3name,b=radius,color=r3color,label=r3label,pos=list(r3))
        r1string='nb_spheres, '+'('+'name'+' '+newr1name+')'
        r2string='nb_spheres, '+'('+'name'+' '+newr2name+')'
        r3string='nb_spheres, '+'('+'name'+' '+r3name+')'
        cmd.show(r1string)
        cmd.show(r2string)
        cmd.show(r3string)
        
        
        
    
    


def GrabMonopoleCharge(atomindex):
    indextotypenum=GenerateIndexToTinkerTypeDic()
    typenum=indextotypenum[atomindex]
    keyfilepath=ReadTinkerFinalKeyFilePath()
    temp=open(keyfilepath,'r')
    results=temp.readlines()
    temp.close()
    dipolevec=[]
    for lineidx in range(len(results)):
        line=results[lineidx]
        if 'multipole' in line:
            linesplit=line.split()
            atomtype=linesplit[1]
            if atomtype==str(typenum):
                monopole=float(linesplit[-1])
    return monopole

def ShowTinkerMonopoleCharge(atomindex):
    monopole=GrabMonopoleCharge(atomindex)
    lab=str(atomindex)+','+str(monopole)
    cmd.select("index "+str(atomindex))
    cmd.label("sele",lab)
    

def ShowTinkerDipoleVector(atomindex):
    indextocoords=GenerateIndexToCoordsDic()
    coords=np.array(indextocoords[atomindex])
    dipolevec=GrabDipoleVector(atomindex)
    cgo_arrow(coords,dipolevec,color='black',lab='Dipole')

def ShowTinkerDipoleCharges(atomindex):
    dipolevec=GrabDipoleVector(atomindex)
    d=np.array(dipolevec)
    indextocoords=GenerateIndexToCoordsDic()
    coords=np.array(indextocoords[atomindex])
    vec=d-coords
    uvec=vec/np.linalg.norm(vec)
    monopole=GrabMonopoleCharge(atomindex)
    disp=.2
    r1=disp*uvec+coords
    r2=-disp*uvec+coords
    r1name=str(atomindex)+'_'+'PositiveDipoleCharge'
    r2name=str(atomindex)+'_'+'NegativeDipoleCharge'
    r1label=str(monopole)
    r2label='-'+str(monopole)
    r1color='red'
    r2color='blue'
    radius=.05
    cmd.pseudoatom(r1name,name=r1name,b=radius,color=r1color,label=r1label,pos=list(r1))
    cmd.pseudoatom(r2name,name=r2name,b=radius,color=r2color,label=r2label,pos=list(r2))
    r1string='nb_spheres, '+'('+'name'+' '+r1name+')'
    r2string='nb_spheres, '+'('+'name'+' '+r2name+')'
    cmd.show(r1string)
    cmd.show(r2string)
    return r1,r2,r1name,r2name,uvec

def GrabBondParameters(classnumber1,classnumber2,prmfilepath):
    bondline=None
    temp=open(prmfilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
         if 'bond' in line:
            linesplit=line.split()
            if len(linesplit)>=3:
                clses=linesplit[1:3]
                if classnumber1 in clses and classnumber2 in clses:
                    bondline=line
                    return bondline
    return bondline

def GrabAngleParameters(classnumber1,classnumber2,classnumber3,prmfilepath):
    angleline=None
    temp=open(prmfilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'angle' in line:
            linesplit=line.split()
            if len(linesplit)>=4:
                clses=linesplit[1:4]
                if classnumber1 in clses and classnumber2 in clses and classnumber3 in clses:
                    angleline=line
                    return angleline
    return angleline

def GrabStretchBendParameters(classnumber1,classnumber2,classnumber3,prmfilepath):
    strbndline=None
    temp=open(prmfilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'strbnd' in line:
            linesplit=line.split()
            if len(linesplit)>=4:
                clses=linesplit[1:4]
                if classnumber1 in clses and classnumber2 in clses and classnumber3 in clses:
                    strbndline=line
                    return strbndline
    return strbndline

def GrabTorsionParameters(classnumber1,classnumber2,classnumber3,classnumber4,prmfilepath):
    torsionline=None
    temp=open(prmfilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'torsion' in line:
            linesplit=line.split()
            if len(linesplit)>=5:
                tors=linesplit[1:5]
                if classnumber1 in tors and classnumber2 in tors and classnumber3 in tors and classnumber4 in tors:
                    torsionline=line
                    return torsionline
    return torsionline

def PrintTinkerTorsionParameters(atomindex1,atomindex2,atomindex3,atomindex4):
    pathtotinkxyz=ReadTinkerFinalXYZFilePath()
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    prmfilepath=GrabParameterFilePathFromTinkerKey()
    typenumber1=indextotypenum[atomindex1]
    typenumber2=indextotypenum[atomindex2]
    typenumber3=indextotypenum[atomindex3]
    typenumber4=indextotypenum[atomindex4]
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    classnumber1=atomtypetoatomclass[typenumber1]
    classnumber2=atomtypetoatomclass[typenumber2]
    classnumber3=atomtypetoatomclass[typenumber3]
    classnumber4=atomtypetoatomclass[typenumber4]
    torsionline=GrabTorsionParameters(classnumber1,classnumber2,classnumber3,classnumber4,keyfilepath)
    if torsionline==None:
         torsionline=GrabTorsionParameters(classnumber1,classnumber2,classnumber3,classnumber4,prmfilepath)
    if torsionline==None:
        torsionline='Parameters Not Found'
    print(torsionline)
    return torsionline
             
def PrintTinkerAngleParameters(atomindex1,atomindex2,atomindex3):
    pathtotinkxyz=ReadTinkerFinalXYZFilePath()
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    prmfilepath=GrabParameterFilePathFromTinkerKey()
    typenumber1=indextotypenum[atomindex1]
    typenumber2=indextotypenum[atomindex2]
    typenumber3=indextotypenum[atomindex3]
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    classnumber1=atomtypetoatomclass[typenumber1]
    classnumber2=atomtypetoatomclass[typenumber2]
    classnumber3=atomtypetoatomclass[typenumber3]
    angleline=GrabAngleParameters(classnumber1,classnumber2,classnumber3,keyfilepath)
    if angleline==None:
         angleline=GrabAngleParameters(classnumber1,classnumber2,classnumber3,prmfilepath)
    if angleline==None:
        angleline='Parameters Not Found'
    print(angleline)
    return angleline

def PrintTinkerStretchBendParameters(atomindex1,atomindex2,atomindex3):
    pathtotinkxyz=ReadTinkerFinalXYZFilePath()
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    prmfilepath=GrabParameterFilePathFromTinkerKey()
    typenumber1=indextotypenum[atomindex1]
    typenumber2=indextotypenum[atomindex2]
    typenumber3=indextotypenum[atomindex3]
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    classnumber1=atomtypetoatomclass[typenumber1]
    classnumber2=atomtypetoatomclass[typenumber2]
    classnumber3=atomtypetoatomclass[typenumber3]
    strbndline=GrabStretchParameters(classnumber1,classnumber2,classnumber3,keyfilepath)
    if strbndline==None:
         strbndline=GrabStretchBendParameters(classnumber1,classnumber2,classnumber3,prmfilepath)
    if strbndline==None:
        strbndline='Parameters Not Found'
    print(strbndline)
    return strbndline
          
def PrintTinkerBondParameters(atomindex1,atomindex2):
    pathtotinkxyz=ReadTinkerFinalXYZFilePath()
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    prmfilepath=GrabParameterFilePathFromTinkerKey()
    typenumber1=indextotypenum[atomindex1]
    typenumber2=indextotypenum[atomindex2]
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    classnumber1=atomtypetoatomclass[typenumber1]
    classnumber2=atomtypetoatomclass[typenumber2]
    bondline=GrabBondParameters(classnumber1,classnumber2,keyfilepath)
    if bondline==None:
         bondline=GrabBondParameters(classnumber1,classnumber2,prmfilepath)
    if bondline==None:
        bondline='Parameters Not Found'
    print(bondline)
    return bondline                   

def PrintTinkerPolarizationParameters(atomindex):
    pathtotinkxyz=ReadTinkerFinalXYZFilePath()
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    typenumber=indextotypenum[atomindex]
    prmfilepath=GrabParameterFilePathFromTinkerKey()
    polline=GrabPolarizationParameters(typenumber,keyfilepath)
    if polline==None:
         polline=GrabPolarizationParameters(typenumber,prmfilepath)
    if polline==None:
        polline='Parameters Not Found'
    print(polline)
    return polline

def LabelTinkerPolarizationParameters(atomindex):
    polline=PrintTinkerPolarizationParameters(atomindex)
    linesplit=polline.split()
    prm1=linesplit[2]
    prm2=linesplit[3]
    cmd.select("index "+atomindex)
    lab=atomindex+','+prm1+','+prm2
    cmd.label("sele",lab)

def PrintTinkerVdwParameters(atomindex):
    pathtotinkxyz=ReadTinkerFinalXYZFilePath()
    keyfilepath=ReadTinkerFinalKeyFilePath()
    indextotypenum=GenerateIndexToTinkerTypeDic()
    prmfilepath=GrabParameterFilePathFromTinkerKey()
    typenumber=indextotypenum[atomindex]
    atomtypetoatomclass=GenerateTinkerTypeToTinkerClassDic()
    classnumber=atomtypetoatomclass[typenumber]
    vdwline=GrabVdwParameters(classnumber,keyfilepath)
    if vdwline==None:
         vdwline=GrabVdwParameters(classnumber,prmfilepath)
    if vdwline==None:
        vdwline='Parameters Not Found'
    print(vdwline)
    return vdwline

def LabelTinkerVdwParameters(atomindex):
    vdwline=PrintTinkerVdwParameters(atomindex)
    linesplit=vdwline.split()
    prm1=linesplit[2]
    prm2=linesplit[3]
    cmd.select("index "+atomindex)
    lab=atomindex+','+prm1+','+prm2
    cmd.label("sele",lab)
cmd.extend("LabelTinkerTypeNumbers",LabelTinkerTypeNumbers)
cmd.extend("LabelIndexAndTinkerTypeNumbers",LabelIndexAndTinkerTypeNumbers)
cmd.extend("LabelIndexAndTinkerTypeAndClassNumbers",LabelIndexAndTinkerTypeAndClassNumbers)
cmd.extend("LabelTinkerClassNumbers",LabelTinkerClassNumbers)
cmd.extend("PrintTinkerPolarizationParameters",PrintTinkerPolarizationParameters)
cmd.extend("PrintTinkerVdwParameters",PrintTinkerVdwParameters)
cmd.extend("PrintTinkerBondParameters",PrintTinkerBondParameters)
cmd.extend("PrintTinkerAngleParameters",PrintTinkerAngleParameters)
cmd.extend("PrintTinkerStretchBendParameters",PrintTinkerStretchBendParameters)
cmd.extend("PrintTinkerTorsionParameters",PrintTinkerTorsionParameters)
cmd.extend("PrintTinkerMultipoleParameters",PrintTinkerMultipoleParameters)
cmd.extend("LabelTinkerPolarizationParameters",LabelTinkerPolarizationParameters)
cmd.extend("LabelTinkerVdwParameters",LabelTinkerVdwParameters)
cmd.extend("DefineTinkerFinalXYZFilePath",DefineTinkerFinalXYZFilePath)
cmd.extend("DefineTinkerXYZFilePath",DefineTinkerXYZFilePath)
cmd.extend("DefineTinkerKeyFilePath",DefineTinkerKeyFilePath)
cmd.extend("DefineTinkerFinalKeyFilePath",DefineTinkerFinalKeyFilePath)
cmd.extend("DefinePoltypeFilePath",DefinePoltypeFilePath)
cmd.extend("ChangeTinkerVdwParameters",ChangeTinkerVdwParameters)
cmd.extend("ChangeTinkerPolarizationParameters",ChangeTinkerPolarizationParameters)
cmd.extend("ChangeTinkerBondParameters",ChangeTinkerBondParameters)
cmd.extend("ChangeTinkerAngleParameters",ChangeTinkerAngleParameters)
cmd.extend("ChangeTinkerStretchBendParameters",ChangeTinkerStretchBendParameters)
cmd.extend("ChangeTinkerTorsionParameters",ChangeTinkerTorsionParameters)
cmd.extend("AddTinkerVdwParameters",AddTinkerVdwParameters)
cmd.extend("AddTinkerPolarizationParameters",AddTinkerPolarizationParameters)
cmd.extend("AddTinkerBondParameters",AddTinkerBondParameters)
cmd.extend("AddTinkerAngleParameters",AddTinkerAngleParameters)
cmd.extend("AddTinkerStretchBendParameters",AddTinkerStretchBendParameters)
cmd.extend("AddTinkerTorsionParameters",AddTinkerTorsionParameters)
cmd.extend("ShowTinkerMultipoleFrame",ShowTinkerMultipoleFrame)
cmd.extend("DefineTinkerBINPath",DefineTinkerBINPath)
cmd.extend("PrintTinkerTotalCharge",PrintTinkerTotalCharge)
cmd.extend("PrintTinkerComponentEnergy",PrintTinkerComponentEnergy)
cmd.extend("AssumePoltypeDefaults",AssumePoltypeDefaults)
cmd.extend("ShowTinkerDipoleCharges",ShowTinkerDipoleCharges)
cmd.extend("ShowTinkerDipoleVector",ShowTinkerDipoleVector)
cmd.extend("ShowTinkerQuadrupoleCharges",ShowTinkerQuadrupoleCharges)
cmd.extend("ShowTinkerMonopoleCharge",ShowTinkerMonopoleCharge)
