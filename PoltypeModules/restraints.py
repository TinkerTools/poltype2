import keyfilemodifications as keymods
import math
import numpy as np
import mdtraj as md
import sys
import os
import time
import scipy.optimize
from scipy.spatial import distance_matrix
from scipy.cluster.hierarchy import fcluster, leaders, linkage
from scipy.spatial.distance import pdist
from collections import Counter, defaultdict
from openbabel import openbabel
import itertools

def ComputeCOM(poltype,atomidxtovecdic,atomidxtomassdic):
    num=np.array([0,0,0],dtype='float64')
    temp=np.array([0,0,0],dtype='float64')
    totmass=0
    for atomidx in atomidxtovecdic.keys():
        mass=atomidxtomassdic[atomidx]
        vec=atomidxtovecdic[atomidx]
        num+=mass*np.array(vec)
        temp+=np.array(vec)
        totmass+=mass
    COM=num/totmass
    return COM


def AddHarmonicRestrainGroupTermsToKeyFile(poltype,keyfilename,teatherdistances,restraintconstant,restrainatomgroup1,restrainatomgroup2):
    count=1
    for i in range(len(restrainatomgroup1)):
        teatherdist=teatherdistances[i]
        first=restrainatomgroup1[i]
        second=restrainatomgroup2[i]
        group1string='group '+str(count)+' '
        group2string='group '+str(count+1)+' '

        for num in first:
            group1string+=str(num)+' '
        keymods.AddKeyWord(poltype,keyfilename,group1string+'\n')       
        for num in second:
            group2string+=str(num)+' '
        keymods.AddKeyWord(poltype,keyfilename,group2string+'\n')
        if poltype.flatbotrest==False:
            restrainstring='restrain-groups '+str(count)+' '+ str(count+1)+' '+str(restraintconstant)+' '+str(teatherdist)+' '+str(teatherdist)
        else:
            restrainstring='restrain-groups '+str(count)+' '+ str(count+1)+' '+str(restraintconstant)+' '+'0'+' '+str(teatherdist)
        keymods.AddKeyWord(poltype,keyfilename,restrainstring+'\n')
        count+=2

    
        
    
def AverageCOMGroups(poltype,filename,grp1,grp2):
    try:
        t = md.load_arc(filename)
    except IOError:
        t = md.load(filename)

    g1COMtog2COMdistarray=[]
    g1atomidxtovecdic={}
    g1atomidxtomassdic={} 
    g2atomidxtovecdic={}
    g2atomidxtomassdic={}
    for framecount in range(t.n_frames):
        for a in t.topology.atoms:
            zeroindex=a.index    
            atomidx=zeroindex+1
            atomsymb=a.element.symbol
            mass=poltype.elementsymtomass[atomsymb]
            vec=t.xyz[framecount,zeroindex,:]*10 # convert nm to Angstrom
            if atomidx in grp1:
                g1atomidxtovecdic[atomidx]=vec
                g1atomidxtomassdic[atomidx]=mass
            elif atomidx in grp2:
                g2atomidxtovecdic[atomidx]=vec
                g2atomidxtomassdic[atomidx]=mass

        if len(g1atomidxtovecdic.keys())==len(grp1) and len(g2atomidxtovecdic.keys())==len(grp2):
            g1COM=ComputeCOM(poltype,g1atomidxtovecdic,g1atomidxtomassdic)
            g2COM=ComputeCOM(poltype,g2atomidxtovecdic,g2atomidxtomassdic)
            dispx=g1COM[0]-g2COM[0]
            dispy=g1COM[1]-g2COM[1]
            dispz=g1COM[2]-g2COM[2]
            dist=np.sqrt(dispx**2+dispy**2+dispz**2)
            g1COMtog2COMdistarray.append(dist)
            atomidxtovecdic={}
            atomidxtomassdic={}
    COM=round(np.mean(g1COMtog2COMdistarray),2)
    if COM<1:
        x=2
    else:
        x=COM+1
    return x



def GroupRestraintFreeEnergyFix(poltype,distances):
    fo=poltype.distancerestraintconstant
    pi=np.pi
    gasconst=8.314/4184
    temp=poltype.equilibriatescheme[-1]
    kt = temp * gasconst
    stdcon = (1.0*10**27) / (6.02214 * 10**23)
    if poltype.flatbotrest==True:
        ri=0
        fi=1
    else:
        ri=poltype.restraintdistance
        fi=poltype.distancerestraintconstant
    poltype.rescorrection=0
    for dist in distances:
         ro=dist
         v1 = 2.0*pi*ri*(-2.0+np.exp(-ri**2*fi/kt))*kt/fi + np.sqrt(kt*(pi/fi)**3)*(2.0*fi*ri*ri+kt)*math.erf(ri*np.sqrt(fi/kt))
         v2 = (4.0*pi/3.0) * (ro**3-ri**3)
         v3 = np.sqrt(kt*(pi/fo)**3) * (2.0*fo*ro*ro+kt+4.0*ro*np.sqrt(kt*fo/pi))
         vol = v1 + v2 + v3
         dv1 = 2.0*pi*ri**3*np.exp(-ri**2*fi/kt)/temp + 2.0*pi*ri*(-2.0+np.exp(-ri**2*fi/kt))*kt/(fi*temp) + 0.5*np.sqrt((pi/fi)**3)*np.sqrt(kt)*(2.0*ri**2*fi+kt)*math.erf(ri*np.sqrt(fi/kt))/temp - pi*ri*np.exp(-ri**2*fi/kt)*(2.0*ri**2*fi+kt)/(fi*temp) + np.sqrt((kt*pi/fi)**3)*math.erf(ri*np.sqrt(fi/kt))/temp
         dv2 = 0.0
         dv3 = np.sqrt(kt*(pi/fo)**3)*fo*ro*ro/temp + 4.0*kt*(pi/fo)*ro/temp + 1.5*np.sqrt((kt*pi/fo)**3)/temp
         dvol = dv1 + dv2 + dv3
         dg = -kt * np.log(vol/stdcon)
         ds = -dg/temp + kt*dvol/vol
         dh = dg + temp*ds
         poltype.rescorrection+=dg
    return dg,dh,ds


def GrabIndexInfo(poltype,fxyz):
    try:
        t = md.load_arc(fxyz)
    except IOError:
        t = md.load(fxyz)
    
    indextomass={}
    indextosym={}
    indextovec={}
    for a in t.topology.atoms:
        zeroindex=a.index    
        atomidx=zeroindex+1
        atomsymb=a.element.symbol
        mass=poltype.elementsymtomass[atomsymb]
        vec=t.xyz[0,zeroindex,:]
        vec=vec*10
        indextovec[zeroindex+1]=vec
        indextomass[zeroindex+1]=mass
        indextosym[zeroindex+1]=atomsymb
    return indextomass,indextovec,indextosym,t


def GrabXYZAndMass(poltype,indextomass,indextovec,indices):
    atomidxtomassdic={}
    atomidxtovecdic={}
    for index in indices:
        atomidxtomassdic[index]=indextomass[index]
        atomidxtovecdic[index]=indextovec[index]


    return atomidxtomassdic,atomidxtovecdic

def FindClosestLigandAtomToCOM(poltype,ligandcom,atomidxtovecdic,indextosym):
    disttoatomidx={}
    for atomidx,vec in atomidxtovecdic.items():
        symb=indextosym[atomidx]
        if symb!='H':
            dist=np.linalg.norm(vec-ligandcom)
            disttoatomidx[dist]=atomidx
    mindist=min(disttoatomidx.keys())
    ligandatomindex=disttoatomidx[mindist]

    return ligandatomindex


def ComputeAverageDistanceFromAtom(poltype,comb,indextovec,ligandatomindex):
    ref=indextovec[ligandatomindex]
    distances=[]
    for index in comb:
        vec=indextovec[index]
        dist=np.linalg.norm(vec-ref)
        distances.append(dist)
    averagedist=np.mean(distances)
    return averagedist


def ComputeIdealGroupRestraints(poltype,ligandindices,fxyz,nmin=4,nmax=8,sele='CA',annihilate=False):
    COMthresh=1
    maxatomcutoff=20
    indextomass,indextovec,indextosym,t=GrabIndexInfo(poltype,fxyz)
    atomidxtomassdic,atomidxtovecdic=GrabXYZAndMass(poltype,indextomass,indextovec,ligandindices)
    ligandcom=ComputeCOM(poltype,atomidxtovecdic,atomidxtomassdic)       
    ligandatomindex=FindClosestLigandAtomToCOM(poltype,ligandcom,atomidxtovecdic,indextosym)
    
    if sele=='protein':
        indices=poltype.proteinindices[:]
    else:
        indices = t.topology.select('name %s'%(sele))
        indices=[i+1 for i in indices]
    newindices=[i for i in indices if i not in ligandindices and i not in poltype.hetatmindices]
    ref=indextovec[ligandatomindex]
    disttoindex={}
    for index in newindices:
        vec=indextovec[index]
        dist=np.linalg.norm(vec-ref)
        disttoindex[dist]=index
    
    sortdists=sorted(disttoindex.keys())
    dists=sortdists[:maxatomcutoff+1]
    allowedindices=[disttoindex[dist] for dist in dists]
    allcombs=[]
    for n in range(nmin,nmax+1):
        combs=itertools.combinations(allowedindices, n)
        for comb in combs:
            allcombs.append(list(comb))
    combtoaveragedistance={}
    for comb in allcombs:
        averagedistance=ComputeAverageDistanceFromAtom(poltype,comb,indextovec,ligandatomindex)
        combtoaveragedistance[tuple(comb)]=averagedistance
    sortedcombtoaveragedistance=dict(sorted(combtoaveragedistance.items(), key=lambda item: item[1]))
    foundres=False
    for comb,avgdist in sortedcombtoaveragedistance.items():
        atomidxtomassdic,atomidxtovecdic=GrabXYZAndMass(poltype,indextomass,indextovec,comb)
        grpcom=ComputeCOM(poltype,atomidxtovecdic,atomidxtomassdic)
        dist=np.linalg.norm(grpcom-ref)
        if dist<COMthresh:
            foundres=True
            poltype.restrainatomgroup2.append(list(comb))
            if annihilate==True:
                poltype.annihilaterestrainatomgroup2.append(list(comb))

            break
        if nmin==1 and len(comb)==1: # last attempt nothing else worked...
            foundres=True
            poltype.restrainatomgroup2.append(list(comb))
            if annihilate==True:
                poltype.annihilaterestrainatomgroup2.append(list(comb))

            break

    if foundres==False:
        raise ValueError('Failed to find ideal host-guest restraint, nmin='+str(nmin)+' , sele='+str(sele))
    else:
        poltype.restrainatomgroup1.append([ligandatomindex])
        if annihilate==True:
            poltype.annihilaterestrainatomgroup1.append([ligandatomindex])
