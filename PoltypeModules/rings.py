from itertools import combinations
import torsiongenerator as torgen
import optimization as opt
import os
import openbabel
import numpy
from itertools import product


def NonAromaticRingAtomicIndices(poltype,mol):
    sssr = mol.GetSSSR()
    atomindices=[]
    for ring in sssr:
        ringatomindices=GrabRingAtomIndices(poltype,mol,ring)
        if ring.IsAromatic()==False:
            atomindices.append(ringatomindices)        

    return atomindices

def GrabRingAtomIndices(poltype,mol,ring):
    ringatomindices=[]
    atomiter=openbabel.OBMolAtomIter(mol)
    for atom in atomiter:
        atomidx=atom.GetIdx()
        if ring.IsInRing(atomidx)==True:
            ringatomindices.append(atomidx)
    return ringatomindices


def NonAromaticRingTorsions(poltype,torsions,atomindices):
    nonarotorsions=[]
    nonarotorsionsflat=[]
    for ring in atomindices:
        ringtors=[]
        for torsion in torsions:
            nonaro=isTorsionInNonAromaticRing(poltype,torsion,ring)
            if nonaro==True:
                ringtors.append(torsion)
                nonarotorsionsflat.append(torsion)
        nonarotorsions.append(ringtors)
    return nonarotorsions,nonarotorsionsflat

def isTorsionInNonAromaticRing(poltype,torsion,ring):
    firstatomindex=torsion[0]
    secondatomindex=torsion[1]
    thirdatomindex=torsion[2]
    fourthatomindex=torsion[3]
    nonaro=False
    if firstatomindex in ring and secondatomindex in ring and thirdatomindex in ring and fourthatomindex in ring:
        nonaro=True
    return nonaro

def TotalParametersToFitForNonAromaticRing(poltype,ringtors): # symmetry class key, not atom indices
    return len(ringtors)*poltype.foldnum+1 #number of torsions x # cosine terms +1 extra profile shift parameter

def TotalDatapointsForNonAromaticRing(poltype,numparameters):
    return numparameters+1 # +1 so dont overfit

def AllPossiblePuckeringLocationsForRing(poltype,ringtors,tortoneighbtors):
    numbertors=len(ringtors)-3 
    combs=list(combinations(ringtors,numbertors))
    finalcombs=[]
    for comb in combs:
        goodcomb=True
        for i in range(len(comb)):
            tor=comb[i]
            neighbtors=tortoneighbtors[tuple(tor)]
            othertors = [x for j,x in enumerate(comb) if j!=i]
            for otor in othertors:
                if otor not in neighbtors:
                    goodcomb=False
        if goodcomb==True:
            finalcombs.append(comb)
    return finalcombs

def NeighboringTorsion(poltype,ringtors,mol):
    tortoneighbtors={}
    for tor in ringtors:
        a,b,c,d=tor[:]
        finaltors=[]
        ntors=IterativeOverNeighborsOfEndTorsionAtoms(poltype,a,ringtors,mol)
        finaltors.extend(ntors)
        ntors=IterativeOverNeighborsOfEndTorsionAtoms(poltype,d,ringtors,mol)
        finaltors.extend(ntors)
        tortoneighbtors[tuple(tor)]=finaltors
    return tortoneighbtors

def IterativeOverNeighborsOfEndTorsionAtoms(poltype,atomidx,ringtors,mol):
    finaltors=[]
    atom=mol.GetAtom(atomidx)        
    atomatomiter=openbabel.OBAtomAtomIter(atom)
    for natom in atomatomiter:
        natomidx=natom.GetIdx()
        ntors=CheckForNeighboringTorsions(poltype,natomidx,ringtors)
        finaltors.extend(ntors)
    return finaltors
    
def CheckForNeighboringTorsions(poltype,natomidx,ringtors):
    ntors=[]
    for tor in ringtors:
        a,b,c,d=tor[:]
        if a==natomidx or b==natomidx or c==natomidx or d==natomidx:
            ntors.append(tor)
    return ntors


def DatapointsPerTorsionForNonArtomaticRing(poltype,ringtors,totaldatapoints):
    return int(totaldatapoints/len(ringtors))

def PhasePerTorsionForNonAromaticRing(poltype,dataptspertorsion,maxrange):
    return int(maxrange/dataptspertorsion)

def UpdateAngleIncrement(poltype,phase,torsion):
    a,b,c,d=torsion[0:4]
    key=str(b)+' '+str(c)
    poltype.rotbndtoanginc[key]=phase

def UpdateMaxRange(poltype,torsion,maxrange):
    a,b,c,d=torsion[0:4]
    key=str(b)+' '+str(c)
    poltype.rotbndtomaxrange[key]=maxrange

def UpdateTorsionSets(poltype,nonarotors):
    for tor in nonarotors:
        torset=[tor]
        if torset in poltype.torlist:
            index=poltype.torlist.index(torset)
            del poltype.torlist[index]
    poltype.torlist.append(nonarotors)


def UpdateVariableTorsions(poltype,nonarotors):
    torset=nonarotors
    poltype.torsettovariabletorlist[tuple(torset)]=torset

def DetermineMaxRanges(poltype,torset,optmol,bondtopology):
    phaseangles=[0]*len(torset)
    if poltype.use_gaus==False and poltype.use_gausoptonly==False:
        prefix='%s-opt-' % (poltype.molecprefix)
        postfix='-opt.xyz' 
        prevstrctfname=torgen.GenerateFilename(poltype,torset,phaseangles,prefix,postfix,optmol)
        cmd = 'cp ../%s %s' % (poltype.logoptfname.replace('.log','.xyz'),prevstrctfname)
        poltype.call_subsystem(cmd,True)

    else:
        prefix='%s-opt-' % (poltype.molecprefix)
        postfix='.log' 
        prevstrctfname=torgen.GenerateFilename(poltype,torset,phaseangles,prefix,postfix,optmol)
        # copy *-opt.log found early by Gaussian to 'prevstrctfname'
        cmd = 'cp ../%s %s' % (poltype.logoptfname,prevstrctfname)
        poltype.call_subsystem(cmd,True)


    variabletorlist=poltype.torsettovariabletorlist[tuple(torset)]
    phaselists=[]
    for tor in torset:
        phaselist=range(0,50,10) # just try fine grid, remove points that dont work 
        phaselists.append(phaselist)
    phaseanglelist=numpy.array(list(product(*phaselists)))
    designatexyz='_determine_maxrange'
    keybase=poltype.key4fname
    keybasepath='../'
    failedgridpoints=[]
    if poltype.skipgridsearch==False:
        for phaseangles in phaseanglelist:
            prevstruct = opt.load_structfile(poltype,prevstrctfname)
            prevstruct = opt.PruneBonds(poltype,prevstruct,bondtopology)
            prevstruct=opt.rebuild_bonds(poltype,prevstruct,optmol)
            prevstrctfname,torxyzfname,newtorxyzfname,keyfname=torgen.tinker_minimize(poltype,torset,optmol,variabletorlist,phaseangles,poltype.torsionrestraint,prevstruct,designatexyz,keybase,keybasepath)
            toralzfname = os.path.splitext(torxyzfname)[0] + '.alz'
            torgen.tinker_analyze(poltype,newtorxyzfname,keyfname,toralzfname)
            term=torgen.AnalyzeTerm(poltype,toralzfname)
            if term==False:
                failedgridpoints.append(phaseangles)
    for angles in failedgridpoints:
        index=phaseanglelist.index(angles)
        del phaseanglelist[index]
    seperate_angle_lists = list(zip(*phaseanglelist)) 
    for i in range(len(seperate_angle_lists)):
        torsion=torset[i]
        angle_list=seperate_angle_lists[i]
        max_angle=max(angle_list)
        min_angle=min(angle_list)
        maxrange=max_angle-min_angle 
        UpdateMaxRange(poltype,torsion,maxrange)

def RefineNonAromaticRingTorsions(poltype,mol,optmol,classkeytotorsionparametersguess):
    if not os.path.isdir('qm-torsion'):
        os.mkdir('qm-torsion')

    os.chdir('qm-torsion')

    bondtopology=torgen.GenerateBondTopology(poltype,optmol)
    atomindices=NonAromaticRingAtomicIndices(poltype,mol)
    nonarotorsions,nonarotorsionsflat=NonAromaticRingTorsions(poltype,poltype.alltorsionslist,atomindices)
    poltype.nonaroringtors=nonarotorsionsflat
    reducednonarotorsions=[]
    for nonarotors in nonarotorsions:
        tortoneighbtors=NeighboringTorsion(poltype,nonarotors,mol)
        combs=AllPossiblePuckeringLocationsForRing(poltype,nonarotors,tortoneighbtors)
        firstcomb=combs[0]
        reducednonarotorsions.append(firstcomb)
        UpdateTorsionSets(poltype,firstcomb)
        UpdateVariableTorsions(poltype,firstcomb)
        torset=tuple(firstcomb)
        DetermineMaxRanges(poltype,torset,optmol,bondtopology)
        numparameters=TotalParametersToFitForNonAromaticRing(poltype,firstcomb)
        datapoints=TotalDatapointsForNonAromaticRing(poltype,numparameters)
        poltype.torsionsettonumptsneeded[torset]=datapoints
        dataptspertorsion=DatapointsPerTorsionForNonArtomaticRing(poltype,firstcomb,datapoints)
        for torsion in torset:
            a,b,c,d=torsion[0:4]
            key=str(b)+' '+str(c)
            maxrange=poltype.rotbndtomaxrange[key]
            phasepertorsion=PhasePerTorsionForNonAromaticRing(poltype,dataptspertorsion,maxrange)
            UpdateAngleIncrement(poltype,phasepertorsion,torsion)
            classkey=torgen.get_class_key(poltype,a,b,c,d)
            classkeysplit=classkey.split()
            cla,clb,clc,cld=classkeysplit[:]
            revclasskey='%s %s %s %s' % (cld, clc, clb, cla)

            if classkey in classkeytotorsionparametersguess.keys(): 
                prms=classkeytotorsionparametersguess[classkey] 
            elif revclasskey in classkeytotorsionparametersguess.keys(): 
                prms=classkeytotorsionparametersguess[revclasskey] 

            poltype.classkeytoinitialprmguess[classkey]=prms

    os.chdir('..')

