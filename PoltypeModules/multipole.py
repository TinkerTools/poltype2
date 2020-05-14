import electrostaticpotential as esp
import time
import os
import sys
import openbabel
import shutil
import re
from collections import deque

def is_in_polargroup(poltype,mol, smarts, bond, f):
    """
    Intent: Check if a given bond is in the group defined by the smarts string
    Input:
        mol: OBMol object
        smarts: Smarts String defining the group (e.g. O-[*])
        bond: bond in question
        f: output file (not used)
    Output:
        True if the bond is in the group, False if not
    Referenced By: gen_peditinfile
    Description: -
    """
    sp = openbabel.OBSmartsPattern()
    openbabel.OBSmartsPattern.Init(sp,smarts)
    sp.Match(mol)
    for i in sp.GetUMapList():
        if ((bond.GetBeginAtomIdx() in i) and
             (bond.GetEndAtomIdx() in i)):
            return True
    return False

def CheckIfAllAtomsSameClass(poltype,classlist):
    allsymm=True
    if len(classlist)>=1:
        firstsymm=poltype.idxtosymclass[classlist[0].GetIdx()]
        for atom in classlist:
            atomidx=atom.GetIdx()
            symmetry=poltype.idxtosymclass[atomidx]
            if symmetry!=firstsymm:
                allsymm=False
    return allsymm

def RemoveFromList(poltype,atomlist,atm):
    newatmlist=[]
    idx=atm.GetIdx()
    for newatm in atomlist:
        if idx!=newatm.GetIdx():
            newatmlist.append(newatm)
    return newatmlist


def AtLeastOneHeavyLf1NeighbNotAtom(poltype,lf1atom,atom):
    foundatleastoneheavy=False
    checkneighbs=[neighb for neighb in openbabel.OBAtomAtomIter(lf1atom)]
    for neighb in checkneighbs:
        if neighb.GetIdx()!=atom.GetIdx() and neighb.GetAtomicNum()!=1:
            foundatleastoneheavy=True
    return foundatleastoneheavy


def AtLeastOneHeavyNeighb(poltype,atom):
    foundatleastoneheavy=False
    checkneighbs=[neighb for neighb in openbabel.OBAtomAtomIter(atom)]
    for neighb in checkneighbs:
        if neighb.GetAtomicNum()!=1:
            foundatleastoneheavy=True
    return foundatleastoneheavy


def GrabHeavyAtomIdx(poltype,lf1atom,atom):
    checkneighbs=[neighb for neighb in openbabel.OBAtomAtomIter(lf1atom)]
    for neighb in checkneighbs:
        if neighb.GetIdx()!=atom.GetIdx() and neighb.GetAtomicNum()!=1:
            return neighb.GetIdx()

def FindUniqueNonRepeatingNeighbors(poltype,nlist):
    d={}
    for b in nlist:
        typenum= poltype.idxtosymclass[b.GetIdx()]
        if typenum not in d.keys():
            d[typenum]=0
        d[typenum]+=1
    uniquetypeneighbsnorepeat=[]
    for typenum in d.keys():
        repeat=d[typenum]
        if repeat==1:
            idx=GrabIndexFromUniqueTypeNumber(poltype,nlist,typenum)
            uniquetypeneighbsnorepeat.append(idx)
    sorteduniquetypeneighbsnorepeat=sorted(uniquetypeneighbsnorepeat,reverse=True)
    return sorteduniquetypeneighbsnorepeat

def GrabIndexFromUniqueTypeNumber(poltype,nlist,typenum):
    for b in nlist:
        idx=b.GetIdx()
        typenumber= poltype.idxtosymclass[idx]
        if typenumber==typenum:
            return idx

def GrabIndexesFromUniqueTypeNumber(poltype,nlist,typenum):
    idxlist=[]
    for b in nlist:
        idx=b.GetIdx()
        typenumber= poltype.idxtosymclass[idx]
        if typenumber==typenum:
            if idx not in idxlist:
                idxlist.append(idx)
    return idxlist




def CheckIfNeighbHasSameType(poltype,a,neighbs):
    check=False
    reftype=poltype.idxtosymclass[a.GetIdx()]
    neighbtypes=[poltype.idxtosymclass[b.GetIdx()] for b in neighbs]
    if reftype in neighbtypes:
        check=True
    return check

def GrabNumberOfConnectedHydrogens(poltype,highestsymneighbnorepeat):
    atomneighbs=[neighb for neighb in openbabel.OBAtomAtomIter(highestsymneighbnorepeat)]
    hydnum=0
    for atom in atomneighbs:
        if atom.GetAtomicNum()==1:
            hydnum+=1
    return hydnum


def gen_peditinfile(poltype,mol):
    """
    Intent: Create a file with local frame definitions for each multipole
    These frame definitions are given as input into tinker's poledit
    Input:
        mol: OBMol object
    Output: *-peditin.txt is created
    Referenced By: main
    Description:
    1. Initialize lfzerox, or local-frame-zero-x-component, array
    2. For each atom a
            a. find the two heaviest atoms that a is bound to
               (heaviest defined by their symmetry class;
               the more atoms an atom is bound to, the higher its symmetry class)
               These two atoms define the local frame for atom a
            b. This information is stored in the localframe1 and localframe2 arrays
    3. If atom a is bound to two atoms of the same symm class (that aren't Hydrogens),
       use these two atoms to define the frame
    4. Find the atoms that are only bound to one atom and define their local frame
       based on the local frame of the one atom they are bound to.
    5. Zero out the x-component of the local frame if there is more than one choice for the
       x-component
       Note: If the atom is only bound to one atom, then the array lfzerox is altered
       If the atom is bound to more than one atom, then the array poltype.localframe2 is altered
       Related to how poledit handles the local frames of atoms with valence 1
    6. Define bisectors; i.e. local frames where both components belong to the same sym class
    7. Write out frames to *-peditin.txt
    8. Write out polarizabilities for certain atom types
    9. Define polarizability groups by 'cutting' certain bonds
       The openbabel method, IsRotor() is used to decide whether a bond is cut or not
    """
    lfzerox = [ False ] * mol.NumAtoms()
    atomindextoremovedipquad={} # for methane need to make dipole and quadupole on the carbon zeroed out, will return this for post proccesing the keyfile after poledit is run
    atomindextoremovedipquadcross={}
    atomtypetospecialtrace={} # for H on CH4 need to make sure Qxx=Qyy=-1/2*Qzz
    idxtobisecthenzbool={}
    idxtobisectidxs={}
    idxtotrisecbool={}
    idxtotrisectidxs={}
    # Assign local frame for atom a based on the symmetry classes of the atoms it is bound to
    atomiter=openbabel.OBMolAtomIter(mol)
    poltype.localframe2 = list(poltype.localframe2)
    for a in atomiter:
        idxtobisecthenzbool[a.GetIdx()]=False
        idxtotrisecbool[a.GetIdx()]=False
    for atom in openbabel.OBMolAtomIter(mol):
        atomidx=atom.GetIdx()
        val=atom.GetValence()
        atomneighbs=[neighb for neighb in openbabel.OBAtomAtomIter(atom)]
        uniqueneighbtypes=list(set([poltype.idxtosymclass[b.GetIdx()] for b in atomneighbs]))
        sorteduniquetypeneighbsnorepeat=FindUniqueNonRepeatingNeighbors(poltype,atomneighbs)
        if len(sorteduniquetypeneighbsnorepeat)!=0:
            highestsymneighbnorepeatidx=sorteduniquetypeneighbsnorepeat[0]
            highestsymneighbnorepeat=mol.GetAtom(highestsymneighbnorepeatidx)
            numhyds=GrabNumberOfConnectedHydrogens(poltype,highestsymneighbnorepeat)
            highestsymneighbnorepeatval=highestsymneighbnorepeat.GetValence()
            neighbsofneighb=[neighb for neighb in openbabel.OBAtomAtomIter(highestsymneighbnorepeat)]
            uniqueneighbtypesofhighestsymneighbnorepeat=list(set([poltype.idxtosymclass[b.GetIdx()] for b in neighbsofneighb]))
            neighbsofneighbwithoutatom=RemoveFromList(poltype,neighbsofneighb,atom)
            uniqueneighbtypesofhighestsymneighbnorepeatwithoutatom=list(set([poltype.idxtosymclass[b.GetIdx()] for b in neighbsofneighbwithoutatom]))

            sorteduniquetypeneighbsnorepeatofneighbsofneighbwithoutatom=FindUniqueNonRepeatingNeighbors(poltype,neighbsofneighbwithoutatom)
        if val==1 and CheckIfAllAtomsSameClass(poltype,[neighb for neighb in openbabel.OBAtomAtomIter(atomneighbs[0])]) and atomneighbs[0].GetValence()==4: # then this is like H in Methane, we want Z-only
            poltype.localframe1[atomidx-1]=sorteduniquetypeneighbsnorepeat[0]
            poltype.localframe2[atomidx - 1] = 0
            lfzerox[atomidx - 1]=True
            atomtypetospecialtrace[atomidx]=True
            atomindextoremovedipquadcross[atomidx]=True
        elif CheckIfAllAtomsSameClass(poltype,atomneighbs) and not AtLeastOneHeavyNeighb(poltype,atom) and val==4: # then this is like carbon in Methane, we want Z-only
            idxlist=GrabIndexesFromUniqueTypeNumber(poltype,atomneighbs,uniqueneighbtypes[0])
            poltype.localframe1[atomidx-1]=idxlist[0]
            poltype.localframe2[atomidx - 1] = 0
            lfzerox[atomidx - 1]=True
            atomindextoremovedipquad[atomidx]=True
        elif atom.GetAtomicNum()==7 and CheckIfAllAtomsSameClass(poltype,atomneighbs) and val==3: # then this is like Ammonia and we can use a trisector here which behaves like Z-only
            idxtotrisecbool[atomidx]=True
            trisectidxs=[atm.GetIdx() for atm in atomneighbs]
            idxtotrisectidxs[atomidx]=trisectidxs
            lfzerox[atomidx - 1]=True # need to zero out the x components just like for z-only case
        elif val==1 and atomneighbs[0].GetValence()==3 and CheckIfAllAtomsSameClass(poltype,[neighb for neighb in openbabel.OBAtomAtomIter(atomneighbs[0])]): # then this is like the H on Ammonia and we can use z-then bisector
            poltype.localframe1[atomidx-1]=sorteduniquetypeneighbsnorepeat[0]
            idxtobisecthenzbool[atomidx]=True
            bisectidxs=[atm.GetIdx() for atm in neighbsofneighbwithoutatom]
            idxtobisectidxs[atomidx]=bisectidxs
        elif (val==2 and CheckIfAllAtomsSameClass(poltype,atomneighbs)) or (val==4 and len(uniqueneighbtypes)==2) and len(sorteduniquetypeneighbsnorepeat)==0:  # then this is like middle propane carbon or oxygen in water
            idxlist=GrabIndexesFromUniqueTypeNumber(poltype,atomneighbs,uniqueneighbtypes[0])
            poltype.localframe1[atomidx-1]=-1*idxlist[0]
            poltype.localframe2[atomidx-1]=-1*idxlist[1]
        elif len(uniqueneighbtypes)==2 and CheckIfNeighbHasSameType(poltype,atom,atomneighbs): # this handles, ethane,ethene...., z-only
            poltype.localframe1[atomidx-1]=sorteduniquetypeneighbsnorepeat[0]
            poltype.localframe2[atomidx - 1] = 0
            lfzerox[atomidx - 1]=True
        elif ((len(uniqueneighbtypes)==2 and val==3) or (len(uniqueneighbtypes)==2 and val==4)) and val==highestsymneighbnorepeatval and len(uniqueneighbtypesofhighestsymneighbnorepeat)==2: # then this is like CH3PO3, we want z-onlytrisector would also work, also handles Analine
            poltype.localframe1[atomidx-1]=sorteduniquetypeneighbsnorepeat[0]
            poltype.localframe2[atomidx - 1] = 0
            lfzerox[atomidx - 1]=True
        elif ((val==1 and (highestsymneighbnorepeatval==3)) or ((val==3) and highestsymneighbnorepeatval==1)) and numhyds<=1 and len(uniqueneighbtypes)<=2 and len(uniqueneighbtypesofhighestsymneighbnorepeat)<=2:
            poltype.localframe1[atomidx-1]=sorteduniquetypeneighbsnorepeat[0]
            poltype.localframe2[atomidx - 1] = 0
            lfzerox[atomidx - 1]=True
        elif (((val==3) and (highestsymneighbnorepeatval==4))) and numhyds<=2 and len(uniqueneighbtypes)<=2 and len(uniqueneighbtypesofhighestsymneighbnorepeat)<=3:
            poltype.localframe1[atomidx-1]=sorteduniquetypeneighbsnorepeat[0]
            poltype.localframe2[atomidx - 1] = 0
            lfzerox[atomidx - 1]=True


        elif ((val==4 and len(uniqueneighbtypes)==2 and highestsymneighbnorepeatval==3) or (val==3 and len(uniqueneighbtypes)==2 and highestsymneighbnorepeatval==4)) and len(uniqueneighbtypesofhighestsymneighbnorepeat)==2 and len(uniqueneighbtypesofhighestsymneighbnorepeatwithoutatom)==1:  # then this is like methyl-amine and we can use the two atoms with same symmetry class to do a z-then-bisector
            idxtobisecthenzbool[atomidx]=True
            bisectidxs=[atm.GetIdx() for atm in neighbsofneighbwithoutatom]
            idxtobisectidxs[atomidx]=bisectidxs
            poltype.localframe1[atomidx - 1] = highestsymneighbnorepeatidx

            # now make sure neighboring atom (lf1) also is using z-then-bisector
        else:
            if len(sorteduniquetypeneighbsnorepeat)==1:
                neighboffirstneighbs=[]
                for n in openbabel.OBAtomAtomIter(highestsymneighbnorepeat):
                    neighboffirstneighbs.append(n)
                newneighbs=RemoveFromList(poltype,neighboffirstneighbs,atom)
                newsorteduniquetypeneighbsnorepeat=FindUniqueNonRepeatingNeighbors(poltype,newneighbs)
                sorteduniquetypeneighbsnorepeat+=newsorteduniquetypeneighbsnorepeat
            poltype.localframe1[atomidx - 1]=sorteduniquetypeneighbsnorepeat[0]
            poltype.localframe2[atomidx - 1]=sorteduniquetypeneighbsnorepeat[1]

    # write out the local frames
    iteratom = openbabel.OBMolAtomIter(mol)
    if not os.path.isfile(poltype.peditinfile):
        f = open (poltype.peditinfile, 'w')
        for a in iteratom:
            if not idxtobisecthenzbool[a.GetIdx()] and not idxtotrisecbool[a.GetIdx()]:
                f.write(str(a.GetIdx()) + " " + str(poltype.localframe1[a.GetIdx() - 1]) + " " + str(poltype.localframe2[a.GetIdx() - 1]) + "\n")
            elif idxtobisecthenzbool[a.GetIdx()] and not idxtotrisecbool[a.GetIdx()]:
                bisectidxs=idxtobisectidxs[a.GetIdx()]
                f.write(str(a.GetIdx()) + " " + str(poltype.localframe1[a.GetIdx() - 1]) + " -" + str(bisectidxs[0])+ " -" + str(bisectidxs[1]) + "\n")
            else:
                trisecidxs=idxtotrisectidxs[a.GetIdx()]
                f.write(str(a.GetIdx()) + " -" + str(trisecidxs[0])+ " -" + str(trisecidxs[1]) + " -" + str(trisecidxs[2])+ "\n")



        f.write("\n")
        f.write('A'+'\n')

        #Find aromatic carbon, halogens, and bonded hydrogens to correct polarizability
        iteratom = openbabel.OBMolAtomIter(mol)
        writesection=False
        lines=[]
        for a in iteratom:
            if (a.GetAtomicNum() == 6 and a.IsAromatic()):
                lines.append(str(a.GetIdx()) + " " + str(1.750) + "\n")
                writesection=True
            elif (a.GetAtomicNum() == 9):
                lines.append(str(a.GetIdx()) + " " + str(0.507) + "\n")
                writesection=True
            elif (a.GetAtomicNum() == 17):
                lines.append(str(a.GetIdx()) + " " + str(2.500) + "\n")
                writesection=True
            elif (a.GetAtomicNum() == 35):
                lines.append(str(a.GetIdx()) + " " + str(3.595) + "\n")
                writesection=True
            elif (a.GetAtomicNum() == 1):
                iteratomatom = openbabel.OBAtomAtomIter(a)
                for b in iteratomatom:
                    if (b.GetAtomicNum() == 6 and b.IsAromatic()):
                        lines.append(str(a.GetIdx()) + " " + str(0.696) + "\n")
                        writesection=True
        if writesection:
            for line in lines:
                f.write(line)

        # Carboxylate ion O-
        sp = openbabel.OBSmartsPattern()
        openbabel.OBSmartsPattern.Init(sp,'[OD1]~C~[OD1]')
        sp.Match(mol)
        for ia in sp.GetMapList():
            f.write(str(ia[0]) + " " + str(0.921) + "\n")



        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
        #Define polarizable groups by cutting bonds
        iterbond = openbabel.OBMolBondIter(mol)
        for b in iterbond:
            if (b.IsRotor()):
                cut_bond = True
                # If in this group, then don't cut bond
                cut_bond = cut_bond and (not is_in_polargroup(poltype,mol,'a[CH2][*]', b,f))
                cut_bond = cut_bond and (not is_in_polargroup(poltype,mol,'[#6]O-[#1]', b,f))
                #cut_bond = cut_bond and (not is_in_polargroup(mol,'[#6][#6]~O', b,f))
                # Formamide RC=O
                cut_bond = cut_bond and (not is_in_polargroup(poltype,mol,'[*][CH]=O', b,f))
                # Amide
                cut_bond = cut_bond and (not is_in_polargroup(poltype,mol,'N(C=O)', b,f))
                cut_bond = cut_bond and (not is_in_polargroup(poltype,mol,'C(C=O)', b,f))
                cut_bond = cut_bond and (not is_in_polargroup(poltype,mol,'aN', b,f))
                cut_bond = cut_bond and (not is_in_polargroup(poltype,mol,'O-[*]', b,f))
                #cut_bond = cut_bond and (not is_in_polargroup(mol,'C[NH2]', b,f))
                if (cut_bond):
                    f.write( str(b.GetBeginAtomIdx()) + " " + str(b.GetEndAtomIdx()) + "\n")

        #f.write('\n')
        f.write('\n')
        f.write("N\n")
        f.write("N\n")


        f.flush()
        os.fsync(f.fileno())

        f.close()
    return lfzerox,atomindextoremovedipquad,atomtypetospecialtrace,atomindextoremovedipquadcross



def post_proc_localframes(poltype,keyfilename, lfzerox,atomindextoremovedipquad,atomindextoremovedipquadcross):
    """
    Intent: This method runs after the tinker tool Poledit has run and created an
    initial *.key file. The local frames for each multipole are "post processed".
    Zeroed out x-components of the local frame are set back to their original values
    If certain multipole values were not zeroed out by poltype this method zeroes them out
    Input:
       keyfilename: string containing file name of *.key file
       lfzerox: array containing a boolean about whether the x-component of the local frame
                for a given atom should be zeroed or not. Filled in method 'gen_peditin'.
                'lfzerox' is true for atoms that are only bound to one other atom (valence = 1)
                that have more than one possible choice for the x-component of their local frame
    Output: *.key file is edited
    Referenced By: main
    Description:
    1. Move the original *.key file to *.keyb
    2. Create new empty file *.key
    3. Iterate through the lines of *.keyb.
        a. If poledit wrote out the local frame with an x-component missing or as 0
        Then rewrite it with the original x-component (lf2) found in gen_peditin
        b. If poledit did not zero out the local frame x-component for an atom but
        lfzerox is true, zero out the necessary multipole components manually
    4.  Write the edited multipole information to *.key
    """
    # mv *.key to *.keyb
    tmpfname = keyfilename + "b"
    shutil.move(keyfilename, tmpfname)
    keyfh = open(tmpfname)
    lines = keyfh.readlines()
    newlines = []
    # open *.key which is currently empty
    tmpfh = open(keyfilename, "w")
    mpolelines = 0
    # iterate over lines in the *.key file
    for ln1 in range(len(lines)):
        # skip over the lines containing multipole values
        if mpolelines > 0:
            mpolelines -= 1
            continue
        elif 'multipole' in lines[ln1]:
            # Check what poledit wrote as localframe2
            tmplst = lines[ln1].split()
            if len(tmplst) == 5:
                (keywd,atmidx,lf1,lf2,chg) = tmplst
            elif len(tmplst) == 4:
                (keywd,atmidx,lf1,chg) = tmplst
                lf2 = '0'
            elif len(tmplst) == 6:
                (keywd,atmidx,lf1,lf2,lf3,chg) = tmplst
            # manually zero out components of the multipole if they were not done by poledit
            if lfzerox[int(atmidx) - 1]:
                tmpmp = list(map(float, lines[ln1+1].split()))
                tmpmp[0] = 0
                lines[ln1+1] = '%46.5f %10.5f %10.5f\n' % tuple(tmpmp)
                tmpmp = list(map(float, lines[ln1+3].split()))
                tmpmp[0] = 0
                lines[ln1+3] = '%46.5f %10.5f\n' % tuple(tmpmp)
                tmpmp = list(map(float, lines[ln1+4].split()))
                tmpmp[0] = 0
                tmpmp[1] = 0
                lines[ln1+4] = '%46.5f %10.5f %10.5f\n' % tuple(tmpmp)

            if int(atmidx) in atomindextoremovedipquad.keys():
                tmpmp = list(map(float, lines[ln1+1].split()))
                tmpmp[0] = 0
                tmpmp[1] = 0
                tmpmp[2] = 0
                lines[ln1+1] = '%46.5f %10.5f %10.5f\n' % tuple(tmpmp)
                tmpmp = list(map(float, lines[ln1+2].split()))
                tmpmp[0] = 0
                lines[ln1+2] = '%46.5f\n' % tuple(tmpmp)
                tmpmp = list(map(float, lines[ln1+3].split()))
                tmpmp[0] = 0
                tmpmp[1] = 0
                lines[ln1+3] = '%46.5f %10.5f\n' % tuple(tmpmp)
                tmpmp = list(map(float, lines[ln1+4].split()))
                tmpmp[0] = 0
                tmpmp[1] = 0
                tmpmp[2] = 0
                lines[ln1+4] = '%46.5f %10.5f %10.5f\n' % tuple(tmpmp)

            if int(atmidx) in atomindextoremovedipquadcross.keys():

                tmpmp = list(map(float, lines[ln1+3].split()))
                tmpmp[0] = 0
                lines[ln1+3] = '%46.5f %10.5f\n' % tuple(tmpmp)
                tmpmp = list(map(float, lines[ln1+4].split()))
                tmpmp[0] = 0
                tmpmp[1] = 0
                lines[ln1+4] = '%46.5f %10.5f %10.5f\n' % tuple(tmpmp)

            if len(tmplst) == 4:
                linesplit=re.split(r'(\s+)', lines[ln1])
                newtmplist=linesplit[:len(linesplit)-4]
                newtmplist.append('    0')
                newtmplist.append(linesplit[-4])
                newtmplist.append(linesplit[-3])
                newtmplist.append(linesplit[-2])
                templine=''.join(newtmplist)
                newlines.extend(templine)
                newlines.extend(lines[ln1+1:ln1+5])
            else:
                newlines.extend(lines[ln1:ln1+5])
            mpolelines = 4
        # append any lines unrelated to multipoles as is
        else:
            newlines.append(lines[ln1])
    # write out the new lines to *.key
    for nline in newlines:
        tmpfh.write(nline)
    tmpfh.close()
    keyfh.close()


def post_process_mpoles(poltype,keyfilename, scalelist):
    """
    Intent: Iterate through multipoles in 'keyfilename' and scale them if necessary
    Calls 'scale_multipoles'
    Input:
        keyfilename: file name for key file
        scalelist: structure containing scaling information. Found in process_types
    Output: new *.key file
    Referenced By: main
    Description:
    1. Move old keyfile to tmpfname
    2. open an empty key file
    3. iterate through tmpfname scaling the multipoles if necessary
    4. output new multipoles to the new key file
    """
    tmpfname = keyfilename + "b"
    shutil.move(keyfilename, tmpfname)
    keyfh = open(tmpfname)
    lines = keyfh.readlines()
    newlines = []
    tmpfh = open(keyfilename, "w")
    mpolelines = 0
    for ln1 in range(len(lines)):
        if mpolelines > 0:
            mpolelines -= 1
            continue
        elif 'multipole' in lines[ln1] and len(lines[ln1].split())==5:
            (keywd,symcls,lf1,lf2,chg) = lines[ln1].split()
            newlines.extend(scale_multipoles(poltype,symcls,lines[ln1:ln1+5],scalelist))
            mpolelines = 4
        elif 'multipole' in lines[ln1] and len(lines[ln1].split())==6:
            (keywd,symcls,lf1,lf2,lf3,chg) = lines[ln1].split()
            newlines.extend(scale_multipoles(poltype,symcls,lines[ln1:ln1+5],scalelist))
            mpolelines = 4
        else:
            newlines.append(lines[ln1])

    for nline in newlines:
        tmpfh.write(nline)
    tmpfh.close()
    keyfh.close()

def process_types(poltype,mol):
    """
    Intent: Set up scalelist array for scaling certain multipole values
For alchol, the quadrupole on O and H should be mannually scaled by 0.6. This only applies to OH that connect to sp3 Carbon. Similarly for NH in amine (that connects to a sp3 C), scale the Q by 0.75 or 75%. See JCC 2011 32(5):967-77.
    """
    scalelist = {}
    multipole_scale_dict = {}

    for atm in openbabel.OBMolAtomIter(mol):
        if poltype.idxtosymclass[atm.GetIdx()] not in scalelist:
            scalelist[poltype.idxtosymclass[atm.GetIdx()]] = []
            scalelist[poltype.idxtosymclass[atm.GetIdx()]].append(None)
            scalelist[poltype.idxtosymclass[atm.GetIdx()]].append(None)
            scalelist[poltype.idxtosymclass[atm.GetIdx()]].append(None)
            multipole_scale_dict = {}

    #multipole_scale_dict['[OH][CX4]'] = [2, 0.6]
    #multipole_scale_dict['[NH2][CX4]'] = [2, 0.75]
    for (sckey, scval) in multipole_scale_dict.items():
        sp = openbabel.OBSmartsPattern()
        openbabel.OBSmartsPattern.Init(sp,sckey)
        match=sp.Match(mol)
        for ia in sp.GetUMapList():
            scalelist[poltype.idxtosymclass[ia[0]]][scval[0]] = scval[1]

    return scalelist

def scale_multipoles(poltype,symmclass, mpolelines,scalelist):
    """
    Intent: Scale multipoles based on value in scalelist
    """
    symmclass = int(symmclass)
    if scalelist[symmclass][2]:
        qp1 = float(mpolelines[2])
        qp1 = map(float, mpolelines[2].split())
        qp1 = [ scalelist[symmclass][2] * x for x in qp1 ]
        mpolelines[2] = '%46.5f\n' % qp1[0]
        qp2 = map(float, mpolelines[3].split())
        qp2 = [ scalelist[symmclass][2] * x for x in qp2 ]
        mpolelines[3] = '%46.5f %10.5f\n' % tuple(qp2)
        qp3 = map(float, mpolelines[4].split())
        qp3 = [ scalelist[symmclass][2] * x for x in qp3 ]
        mpolelines[4] = '%46.5f %10.5f %10.5f\n' % tuple(qp3)
    return mpolelines

def rm_esp_terms_keyfile(poltype,keyfilename):
    """
    Intent: Remove unnecessary terms from the key file
    """
    tmpfname = keyfilename + "_tmp"
    cmdstr = "sed -e \'/\(.* none$\|^#\|^fix\|^potential-offset\)/d\' " +  keyfilename + " > " + tmpfname
    poltype.call_subsystem(cmdstr,True)
    shutil.move(tmpfname, keyfilename)
    # Removes redundant multipole definitions created after potential fitting
    keyfh = open(keyfilename)
    lines = keyfh.readlines()
    newlines = []
    tmpfh = open(tmpfname, "w")
    mpolelines = 0
    for ln1 in range(len(lines)):
        if mpolelines > 0:
            mpolelines -= 1
            continue
        if 'multipole' in lines[ln1]:
            nllen = len(newlines)
            for ln2 in range(len(newlines)):
                if len(lines[ln1].split()) <= 1 or len(lines[ln2].split()) <= 1:
                    continue
                if lines[ln1].split()[0] in newlines[ln2].split()[0] and \
                lines[ln1].split()[0] in 'multipole' and \
                lines[ln1].split()[1] in newlines[ln2].split()[1]:
                    newlines[ln2] = lines[ln1]
                    newlines[ln2+1] = lines[ln1+1]
                    newlines[ln2+2] = lines[ln1+2]
                    newlines[ln2+3] = lines[ln1+3]
                    newlines[ln2+4] = lines[ln1+4]
                    mpolelines = 4
                    break
            if mpolelines == 0:
                newlines.append(lines[ln1])
        else:
            newlines.append(lines[ln1])

    for nline in newlines:
        tmpfh.write(nline)
    tmpfh.close()
    keyfh.close()
    shutil.move(tmpfname, keyfilename)

def prepend_keyfile(poltype,keyfilename,optmol,dipole=False):
    """
    Intent: Adds a header to the key file given by 'keyfilename'
    """
    while not os.path.isfile(keyfilename):
        time.sleep(5)
        poltype.WriteToLog('Waiting for '+keyfilename)
    tmpfname = keyfilename + "_tmp"
    tmpfh = open(tmpfname, "w")
    keyfh = open(keyfilename, "r")
    tmpfh.write("parameters " + poltype.paramhead + "\n")

    tmpfh.write("bondterm none\n")
    tmpfh.write("angleterm none\n")
    tmpfh.write("torsionterm none\n")
    tmpfh.write("vdwterm none\n")
    tmpfh.write("fix-monopole\n")
    tmpfh.write("potential-offset 1.0\n\n")
    if poltype.use_gaus:
        logname=poltype.logespfname
    else:
        logname=poltype.logespfname.replace('.log','_psi4.log')
    if dipole:
        qmdipole=esp.GrabQMDipoles(poltype,optmol,logname)
        tmpfh.write('TARGET-DIPOLE'+' '+str(qmdipole[0])+' '+str(qmdipole[1])+' '+str(qmdipole[2])+'\n')

    for line in keyfh:
        tmpfh.write(line)
    shutil.move(tmpfname, keyfilename)


def rotate_list(poltype,l1):
    deq = deque(l1)
    deq.rotate(-1)
    return list(deq)

def gen_gdmain(poltype,gdmainfname,molecprefix,fname,dmamethod):
    """
    Intent: Generate GDMA input file for the molecule
    Input:
        gdmainfname: file name for the gdma input file
        molecprefix: name of the molecule
        fname: file name for the *-dma.fchk file
    Output: *.gdmain file is created
    Referenced By: run_gdma
    Description:
        1. create pointer file dma.fchk to *-dma.fchk
        2. create *.gdmain file and write in all the necessary information for the gdma run
    """
    fnamesym = "dma.fchk"
    try:
        os.symlink(fname,fnamesym)
    except Exception as e:
        print(e)



    #punfname = os.path.splitext(fname)[0] + ".punch"
    punfname = "dma.punch"

    try:
        tmpfh = open(gdmainfname, "w")
    except IOError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
        sys.exit(e.errno)

    tmpfh.write("Title " + molecprefix + " gdmain\n")
    tmpfh.write("\n")
    if poltype.dmamethod=='MP2':
        densitystring='MP2'
    else:
        densitystring='SCF'
    if not poltype.use_gaus or poltype.use_gausoptonly:
        if poltype.dmamethod=='MP2':
            densitystring='CC' # for some reason fchk outputs CC for MP2 density
        tmpfh.write("File " + fnamesym  + " density %s\n"%(densitystring))
    else:
        tmpfh.write("File " + fnamesym  + " density %s\n"%(densitystring))
    tmpfh.write("Angstrom\n")
    tmpfh.write("AU\n")
    tmpfh.write("Multipoles\n")
    tmpfh.write("Switch 0\n") # v1.3, comment out for v2.2
    tmpfh.write("Limit 2\n")
    tmpfh.write("Punch " + punfname + "\n")
    tmpfh.write("Radius H 0.60\n")
    tmpfh.write("Radius S 0.85\n")
    tmpfh.write("Radius P 0.95\n")
    tmpfh.write("Radius Cl 1.10\n")
    tmpfh.write("\n")
    tmpfh.write("Start\n")
    tmpfh.write("\n")
    tmpfh.write("Finish\n")
    tmpfh.close()

def run_gdma(poltype):
    """
    Intent: Runs GDMA to find multipole information
    The GDMA program carries out distributed multipole analysis of the wavefunctions
    calculated by Gaussian (using the *-dma.fchk file)
    Input:
    Output: *.gdmaout is created; this is used as input by poledit, a tinker program
    Referenced By: main
    Description:
    1. Generates the gdma input file by calling 'gen_gdmain'
    2. Runs the following command: gdma < *.gdmain > *.gdmaout
    """

    poltype.WriteToLog("NEED DMA: Executing GDMA")

    if not os.path.isfile(poltype.fckdmafname):
        poltype.fckdmafname = os.path.splitext(poltype.fckdmafname)[0]

    assert os.path.isfile(poltype.fckdmafname), "Error: " + poltype.fckdmafname + " does not exist."+' '+os.getcwd()
    poltype.gdmainfname = poltype.assign_filenames ( "gdmainfname" , ".gdmain")
    gen_gdmain(poltype,poltype.gdmainfname,poltype.molecprefix,poltype.fckdmafname,poltype.dmamethod)

    cmdstr = poltype.gdmaexe + " < " + poltype.gdmainfname + " > " + poltype.gdmafname
    poltype.call_subsystem(cmdstr,True)

    assert os.path.getsize(poltype.gdmafname) > 0, "Error: " + os.getcwd() +' '+os.path.basename(poltype.gdmaexe) + " cannot create .gdmaout file."

def AverageMultipoles(poltype,optmol):
    # gen input file
    gen_avgmpole_groups_file(poltype)
    # call avgmpoles.pl
    avgmpolecmdstr = poltype.avgmpolesexe + " " + poltype.keyfname + " " + poltype.xyzfname + " " + poltype.grpfname + " " + poltype.key2fname + " " + poltype.xyzoutfile + " " + str(poltype.prmstartidx)
    poltype.call_subsystem(avgmpolecmdstr,True)
    prepend_keyfile(poltype,poltype.key2fname,optmol,True)

def gen_avgmpole_groups_file(poltype):
    """
    Intent: Print out *-groups.txt which is a map from symm class to idx
    Also, symm class labels are altered from 1, 2, 3, ... to 401, 402, ...
    (or some other labeling system dependent on 'prmstartidx')
    Input:
    Output:
        *-groups.txt: map from symm group to atom idx
    Referenced By: main
    Description: -
    """
    symclasstoidxlist={}
    for idx,symclass in poltype.idxtosymclass.items():
        if symclass not in symclasstoidxlist.keys():
            symclasstoidxlist[symclass]=[]
        symclasstoidxlist[symclass].append(idx)
    f = open(poltype.grpfname,"w")
    for symclass,idxlist in symclasstoidxlist.items():
        string=str(symclass)+' '
        for idx in idxlist:
            string+=str(idx)+' '
        string+='\n'
        f.write(string)
    f.close()


